// Tables HaplotypeOverviewRows, HaplotypeOverviewBins and
// HaplotypeOverviewClasses: for each reference sample with anchors, a binned
// summary of every haplotype against each of its reference paths, at zoom
// levels a factor of four apart, so a browser draws a chromosome-wide view
// from a few reads.
//
// Each path is walked forward once per reference sample. Its visits to the
// sample's reference nodes place it on a reference path; the stretch between
// two consecutive such visits is an excursion, which replaces the reference
// bp between the two nodes with the path's own bp. The path covers the
// reference from the first of the two nodes to the second, and an excursion
// of `sv_bp` or more in either measure marks the bins it spans as variant. A
// contig reads the reference in either direction, so two visits in the same
// orientation that step on in that direction are an alignment; a turn, a step
// back, a jump to another reference path or one of more than MAX_JUMP bp
// marks both nodes' bins as variant instead.
//
// A haplotype is the set of paths with one sample and phase. Per bin it gets
// one of four classes: absent when none of its paths covers the bin,
// variant when one marked it, partial when they cover less than nine tenths
// of it, and reference otherwise; a variant cell also carries a bucket of how
// many marks the bin holds, which tells haplotypes apart at the coarse levels
// where nearly every bin holds one. Per bin the tables also count the
// haplotypes of each class and the excursions that start in it.

use super::{fit_u32, Anchors, PathSource};
use crate::encode::push_varint;
use crate::strays::{reference_samples, Runner};
use gbz::{support, Orientation, ENDMARKER};

use std::collections::{BTreeMap, HashMap};
use std::sync::Mutex;
use std::time::Instant;

pub const FORMAT: &str = "1";
pub const ABSENT: u8 = 0;
pub const REFERENCE: u8 = 1;
pub const PARTIAL: u8 = 2;
pub const VARIANT: u8 = 3;
const MAX_JUMP: usize = 10_000_000;
const ZOOM: usize = 4;

#[derive(Clone, Copy)]
pub struct Options {
    pub bin: usize,
    pub chunk: usize,
    pub sv_bp: usize,
    pub part_bytes: usize,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord, Debug)]
pub struct Row {
    pub sample: String,
    pub haplotype: usize,
}

pub struct Part {
    pub reference_handle: u32,
    pub level: u32,
    pub chunk: u32,
    pub part: u32,
    pub bytes: Vec<u8>,
}

pub struct Output {
    pub rows: Vec<Row>,
    pub bins: Vec<Part>,
    pub classes: Vec<Part>,
    pub levels: usize,
    pub samples: Vec<String>,
}

impl Output {
    pub fn empty() -> Self {
        Output {
            rows: Vec::new(),
            bins: Vec::new(),
            classes: Vec::new(),
            levels: 0,
            samples: Vec::new(),
        }
    }
}

// The reference nodes of one sample: per node id, the reference path, the
// offset where the node starts and whether the reference visits it reversed.
struct RefMap {
    slots: Vec<u64>,
    handles: Vec<usize>,
    lengths: Vec<usize>,
}

impl RefMap {
    fn build(source: &dyn PathSource, handles: Vec<usize>) -> Self {
        let mut slots = vec![0u64; source.max_node_id() + 1];
        let mut lengths = Vec::new();
        for (i, &handle) in handles.iter().enumerate() {
            let mut pos = source.start(handle, Orientation::Forward);
            let mut offset = 0usize;
            while let Some(current) = pos {
                if current.node == ENDMARKER {
                    break;
                }
                let (len, next) = source.step(current);
                let slot = &mut slots[support::node_id(current.node)];
                if *slot == 0 {
                    *slot = ((i as u64 + 1) << 33)
                        | ((current.node as u64 & 1) << 32)
                        | fit_u32("reference offset", offset) as u64;
                }
                offset += len;
                pos = next;
            }
            lengths.push(offset);
        }
        RefMap {
            slots,
            handles,
            lengths,
        }
    }

    fn get(&self, node_handle: usize) -> Option<(usize, usize, bool)> {
        let slot = self.slots[support::node_id(node_handle)];
        if slot == 0 {
            return None;
        }
        let reversed = (slot >> 32) & 1 != (node_handle as u64 & 1);
        Some((
            (slot >> 33) as usize - 1,
            (slot & 0xffff_ffff) as usize,
            reversed,
        ))
    }
}

// The bins of every reference path of the sample, numbered in one sequence.
struct Grid {
    bin: usize,
    first: Vec<usize>,
    lengths: Vec<usize>,
    total: usize,
}

impl Grid {
    fn new(bin: usize, lengths: &[usize]) -> Self {
        let mut first = Vec::new();
        let mut total = 0;
        for &length in lengths {
            first.push(total);
            total += length.div_ceil(bin).max(1);
        }
        Grid {
            bin,
            first,
            lengths: lengths.to_vec(),
            total,
        }
    }

    fn count(&self, reference: usize) -> usize {
        self.lengths[reference].div_ceil(self.bin).max(1)
    }

    fn at(&self, reference: usize, offset: usize) -> usize {
        self.first[reference] + (offset / self.bin).min(self.count(reference) - 1)
    }

    // The bp of the reference path inside its bin b, less than `bin` at the end.
    fn width(&self, reference: usize, b: usize) -> usize {
        let start = (b - self.first[reference]) * self.bin;
        self.lengths[reference].saturating_sub(start).min(self.bin)
    }
}

#[derive(Default)]
struct Marks {
    // reference, from, to: covered stretches and variant stretches
    spans: Vec<(usize, usize, usize)>,
    variant: Vec<(usize, usize, usize)>,
    // reference, offset, bp: excursions
    sites: Vec<(usize, usize, usize)>,
}

impl Marks {
    fn cover(&mut self, reference: usize, from: usize, to: usize) {
        if let Some(last) = self.spans.last_mut() {
            if last.0 == reference && from <= last.2 && to >= last.1 {
                last.1 = last.1.min(from);
                last.2 = last.2.max(to);
                return;
            }
        }
        self.spans.push((reference, from, to));
    }
}

fn walk(source: &dyn PathSource, refs: &RefMap, path_handle: usize, sv_bp: usize) -> Marks {
    let mut marks = Marks::default();
    // reference, offset, end and orientation of the last reference visit
    let mut last: Option<(usize, usize, usize, bool)> = None;
    let mut off_bp = 0usize;
    let mut pos = source.start(path_handle, Orientation::Forward);
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (len, next) = source.step(current);
        match refs.get(current.node) {
            None => off_bp += len,
            Some((reference, offset, reversed)) => {
                let end = offset + len;
                if let Some((r1, o1, e1, rev1)) = last {
                    // The reference bp between the visits when the path reads
                    // on along the reference in its orientation.
                    let gap = if r1 != reference || rev1 != reversed {
                        None
                    } else if !reversed && offset >= e1 {
                        Some(offset - e1)
                    } else if reversed && end <= o1 {
                        Some(o1 - end)
                    } else {
                        None
                    };
                    match gap {
                        Some(gap) if gap <= MAX_JUMP => {
                            if gap > 0 || off_bp > 0 {
                                let bp = gap.max(off_bp);
                                let (from, to) = if reversed {
                                    (end, o1.max(end + 1))
                                } else {
                                    (e1, offset.max(e1 + 1))
                                };
                                marks.sites.push((reference, from, bp));
                                if bp >= sv_bp {
                                    marks.variant.push((reference, from, to));
                                }
                            }
                            if reversed {
                                marks.cover(reference, offset, e1);
                            } else {
                                marks.cover(reference, o1, end);
                            }
                        }
                        _ => {
                            marks.variant.push((r1, o1, e1));
                            marks.variant.push((reference, offset, end));
                            marks.cover(reference, offset, end);
                        }
                    }
                } else {
                    marks.cover(reference, offset, end);
                }
                last = Some((reference, offset, end, reversed));
                off_bp = 0;
            }
        }
        pos = next;
    }
    marks
}

// Per haplotype, the bp covered in each bin and the variant marks on it.
struct Haplotype {
    covered: Vec<u16>,
    marks: Vec<u8>,
}

// A cell of HaplotypeOverviewClasses: the class in the low two bits and, for
// a variant cell, a bucket of the bin's marks in the high two: 1, 2 to 3, 4
// to 15, and 16 or more.
pub fn cell(class: u8, marks: usize) -> u8 {
    let bucket = match marks {
        0 | 1 => 0,
        2..=3 => 1,
        4..=15 => 2,
        _ => 3,
    };
    class | if class == VARIANT { bucket << 2 } else { 0 }
}

#[derive(Clone, Copy, Default)]
struct Site {
    small: u32,
    sv: u32,
    max_bp: u32,
}

fn add_span(grid: &Grid, covered: &mut [u16], reference: usize, from: usize, to: usize) {
    let mut at = from;
    while at < to {
        let b = grid.at(reference, at);
        let bin_end = ((at / grid.bin) + 1) * grid.bin;
        let stop = to.min(bin_end);
        covered[b] = covered[b].saturating_add((stop - at) as u16);
        at = stop;
    }
}

struct Level {
    bin: usize,
    counts: Vec<usize>,
}

fn class_of(covered: usize, width: usize, marks: usize) -> u8 {
    if covered == 0 {
        ABSENT
    } else if marks > 0 {
        VARIANT
    } else if covered * 10 < width * 9 {
        PARTIAL
    } else {
        REFERENCE
    }
}

fn cut(bytes: Vec<u8>, limit: usize) -> Vec<Vec<u8>> {
    if bytes.is_empty() {
        return vec![bytes];
    }
    bytes.chunks(limit.max(1)).map(|c| c.to_vec()).collect()
}

// `lengths` gives every path's length, so that the levels are the same for
// every reference sample: as many as it takes for the longest anchored path
// to fit one chunk.
pub fn overview(
    runner: &dyn Runner,
    anchors: &Anchors,
    lengths: &[(usize, usize)],
    options: Options,
) -> Output {
    let mut out = Output::empty();
    if options.bin == 0 {
        return out;
    }
    let source = runner.source();
    let paths = source.path_count();
    let longest = lengths
        .iter()
        .filter(|(handle, _)| anchors.reference_paths.contains(handle))
        .map(|&(_, length)| length)
        .max()
        .unwrap_or(0);
    let mut level_bins = vec![options.bin];
    while longest.div_ceil(*level_bins.last().unwrap()) > options.chunk {
        level_bins.push(level_bins.last().unwrap() * ZOOM);
    }
    out.levels = level_bins.len();
    let mut by_name: BTreeMap<Row, Vec<usize>> = BTreeMap::new();
    for handle in 0..paths {
        by_name
            .entry(Row {
                sample: source.sample_of(handle),
                haplotype: source.haplotype_of(handle),
            })
            .or_default()
            .push(handle);
    }
    let mut row_of = vec![0usize; paths];
    for (row, handles) in by_name.values().enumerate() {
        for &handle in handles {
            row_of[handle] = row;
        }
    }
    out.rows = by_name.keys().cloned().collect();
    let haplotypes = out.rows.len();
    for (name, reference) in reference_samples(source, anchors) {
        let started = Instant::now();
        let refs = RefMap::build(source, reference.iter().copied().collect());
        let grid = Grid::new(options.bin, &refs.lengths);
        let accumulators: Vec<Mutex<Haplotype>> = (0..haplotypes)
            .map(|_| {
                Mutex::new(Haplotype {
                    covered: vec![0; grid.total],
                    marks: vec![0; grid.total],
                })
            })
            .collect();
        let sites = Mutex::new(vec![Site::default(); grid.total]);
        runner.run(paths, &|source, handle| {
            let marks = walk(source, &refs, handle, options.sv_bp);
            let mut local: HashMap<usize, Site> = HashMap::new();
            for &(reference, offset, bp) in &marks.sites {
                let site = local.entry(grid.at(reference, offset)).or_default();
                if bp >= options.sv_bp {
                    site.sv += 1;
                } else {
                    site.small += 1;
                }
                site.max_bp = site.max_bp.max(bp.min(u32::MAX as usize) as u32);
            }
            {
                let mut acc = accumulators[row_of[handle]].lock().unwrap();
                for &(reference, from, to) in &marks.spans {
                    add_span(&grid, &mut acc.covered, reference, from, to);
                }
                for &(reference, from, to) in &marks.variant {
                    for b in grid.at(reference, from)
                        ..=grid.at(reference, to.saturating_sub(1).max(from))
                    {
                        acc.marks[b] = acc.marks[b].saturating_add(1);
                    }
                }
            }
            if !local.is_empty() {
                let mut shared = sites.lock().unwrap();
                for (b, site) in local {
                    let s = &mut shared[b];
                    s.small += site.small;
                    s.sv += site.sv;
                    s.max_bp = s.max_bp.max(site.max_bp);
                }
            }
        });
        let accumulators: Vec<Haplotype> = accumulators
            .into_iter()
            .map(|m| m.into_inner().unwrap())
            .collect();
        let sites = sites.into_inner().unwrap();
        let levels: Vec<Level> = level_bins
            .iter()
            .map(|&bin| Level {
                bin,
                counts: refs
                    .lengths
                    .iter()
                    .map(|&l| l.div_ceil(bin).max(1))
                    .collect(),
            })
            .collect();
        let bytes_per_bin = haplotypes.div_ceil(2);
        for (level_index, level) in levels.iter().enumerate() {
            let factor = level.bin / options.bin;
            for (reference, &handle) in refs.handles.iter().enumerate() {
                let mut classes = vec![0u8; level.counts[reference] * bytes_per_bin];
                let mut summaries: Vec<[usize; 7]> = vec![[0; 7]; level.counts[reference]];
                for b in 0..level.counts[reference] {
                    let fine_from = grid.first[reference] + b * factor;
                    let fine_to =
                        (fine_from + factor).min(grid.first[reference] + grid.count(reference));
                    let width: usize = (fine_from..fine_to).map(|f| grid.width(reference, f)).sum();
                    let summary = &mut summaries[b];
                    for f in fine_from..fine_to {
                        summary[4] += (sites[f].small + sites[f].sv) as usize;
                        summary[5] += sites[f].sv as usize;
                        summary[6] = summary[6].max(sites[f].max_bp as usize);
                    }
                    for (row, acc) in accumulators.iter().enumerate() {
                        let covered: usize =
                            (fine_from..fine_to).map(|f| acc.covered[f] as usize).sum();
                        let marks: usize =
                            (fine_from..fine_to).map(|f| acc.marks[f] as usize).sum();
                        let class = class_of(covered, width, marks);
                        summary[class as usize] += 1;
                        classes[b * bytes_per_bin + row / 2] |=
                            cell(class, marks) << (4 * (row % 2));
                    }
                }
                for chunk in 0..level.counts[reference].div_ceil(options.chunk) {
                    let from = chunk * options.chunk;
                    let to = (from + options.chunk).min(level.counts[reference]);
                    let mut bins = Vec::new();
                    for summary in &summaries[from..to] {
                        for &value in summary {
                            push_varint(&mut bins, value);
                        }
                    }
                    for (part, bytes) in cut(bins, options.part_bytes).into_iter().enumerate() {
                        out.bins.push(Part {
                            reference_handle: handle as u32,
                            level: level_index as u32,
                            chunk: fit_u32("chunk", chunk),
                            part: part as u32,
                            bytes,
                        });
                    }
                    let slice = classes[from * bytes_per_bin..to * bytes_per_bin].to_vec();
                    for (part, bytes) in cut(slice, options.part_bytes).into_iter().enumerate() {
                        out.classes.push(Part {
                            reference_handle: handle as u32,
                            level: level_index as u32,
                            chunk: fit_u32("chunk", chunk),
                            part: part as u32,
                            bytes,
                        });
                    }
                }
            }
        }
        eprintln!(
            "Overview of {} haplotypes over {} reference paths of {} in {} levels from {} bp bins, {} class parts, in {:.0} s",
            haplotypes,
            refs.handles.len(),
            name,
            levels.len(),
            options.bin,
            out.classes.len(),
            started.elapsed().as_secs_f64()
        );
        out.samples.push(name);
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn classes_follow_coverage_and_marks() {
        assert_eq!(class_of(0, 100, 0), ABSENT);
        assert_eq!(class_of(0, 100, 1), ABSENT);
        assert_eq!(class_of(100, 100, 0), REFERENCE);
        assert_eq!(class_of(90, 100, 0), REFERENCE);
        assert_eq!(class_of(89, 100, 0), PARTIAL);
        assert_eq!(class_of(89, 100, 1), VARIANT);
        assert_eq!(cell(PARTIAL, 0), PARTIAL);
        assert_eq!(cell(VARIANT, 1), VARIANT);
        assert_eq!(cell(VARIANT, 3), VARIANT | 1 << 2);
        assert_eq!(cell(VARIANT, 15), VARIANT | 2 << 2);
        assert_eq!(cell(VARIANT, 16), VARIANT | 3 << 2);
    }

    #[test]
    fn spans_add_their_bp_to_each_bin_they_cross() {
        let grid = Grid::new(100, &[250, 40]);
        assert_eq!(grid.total, 4);
        assert_eq!(grid.width(0, 2), 50);
        assert_eq!(grid.width(1, 3), 40);
        let mut covered = vec![0u16; grid.total];
        add_span(&grid, &mut covered, 0, 90, 230);
        assert_eq!(covered, vec![10, 100, 30, 0]);
        add_span(&grid, &mut covered, 1, 0, 40);
        assert_eq!(covered[3], 40);
    }

    #[test]
    fn consecutive_covers_merge_in_either_direction() {
        let mut marks = Marks::default();
        marks.cover(0, 0, 10);
        marks.cover(0, 10, 20);
        marks.cover(0, 30, 40);
        marks.cover(1, 0, 5);
        assert_eq!(marks.spans, vec![(0, 0, 20), (0, 30, 40), (1, 0, 5)]);
        let mut reversed = Marks::default();
        reversed.cover(0, 90, 100);
        reversed.cover(0, 80, 90);
        reversed.cover(0, 60, 70);
        assert_eq!(reversed.spans, vec![(0, 80, 100), (0, 60, 70)]);
    }
}
