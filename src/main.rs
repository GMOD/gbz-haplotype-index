use gbz::bwt::BWT;
use gbz::support;
use gbz::{Orientation, Pos, ENDMARKER, GBWT, GBZ};
use gbz_base::{GBZBase, GraphInterface};
use rusqlite::{params, Connection, OpenFlags};
use simple_sds::serialize;

mod encode;
mod overview;
mod strays;

use std::collections::{BTreeMap, BTreeSet};
use std::env;
use std::process;
use std::thread;
use std::time::Instant;

const USAGE: &str = "Usage: gbz-haplotype-index [options] graph.gbz [graph.gbz.db] index.db
       gbz-haplotype-index [options] --from-db graph.gbz.db index.db

Writes a haplotype index for a gbz-base database into a separate file,
index.db. Give graph.gbz.db alongside graph.gbz to check that the two match.

Walks every path in both orientations and writes a sample every --interval bp
(default 4096) into table HaplotypeSamples, plus the path lengths into
HaplotypeLengths. The path start and end are always sampled. The reference
paths, or those of --anchor-sample, take a sample every --reference-interval bp
instead, so a query whose subgraph reaches the reference far from its window
finds a reference sample there.

Anchors: along every reference path (the paths gbz-base indexes for random
access), or only those of one sample with --anchor-sample, one anchor node is
chosen per multiple k of --anchor-spacing: the path's first node for k = 0, and
for k >= 1 the node with the most GBWT positions among those overlapping the
half spacing before k * spacing, the one nearest k * spacing on a tie, so it is
one most haplotypes of the region visit. Table HaplotypeAnchors records each
choice, and every path visit through an anchor node is sampled, in both
orientations. A reader can then list every haplotype passing a reference
position, with its own coordinate, from the rows at one node.
--anchor-spacing 0 writes none.

The index records the graph's path and node counts so the reader catches a
mismatch at open.

Strays: for each reference sample with anchors, table HaplotypeBins lists
the nodes within --stray-context bp of each --stray-bin bp of a reference
path, and beside them the visits to those nodes that the keep route's walks
from anchor visits can miss, as rows the reader walks for the paths it keeps.
A query checks its subgraph against the node lists, so the table holds for a
context of at most --stray-context. Give graph.gbz.db, which holds the
top-level snarls, for the rows to cover a query that fills snarls.
--stray-context 0 writes none.

Overview: for each reference sample with anchors, tables HaplotypeOverviewRows,
HaplotypeOverviewBins and HaplotypeOverviewClasses summarize every haplotype
against each reference path in bins of --overview-bin bp and coarser levels a
factor of four apart: per bin and haplotype one of absent, reference, partial
and variant, and per bin the counts of each class and of the excursions from
the reference, where an excursion of --overview-sv bp or more is a variant.
A browser draws a whole chromosome from these in a few reads.
--overview-bin 0 writes none.

With --from-db the walk reads node records from the database itself, so the
GBZ is not needed. Walking a GBZ uses --threads (default: all cores).

Options:
  --interval BP        bp between samples along a path (default 4096)
  --anchor-spacing BP  bp between anchors along a reference path (default 32768)
  --anchor-sample NAME anchor only the reference paths of this sample
  --reference-interval BP
                       bp between samples along a reference path
                       (default: --interval)
  --stray-context BP   the largest query context the stray rows cover
                       (default 1000)
  --stray-bin BP       bp of a reference path per node list (default 16384)
  --stray-bound BP     bp a walk runs past an anchor visit (default 32768)
  --stray-gap BP       bp between strays that starts a new row (default 1024)
  --overview-bin BP    bp of a reference path per overview bin (default 4096)
  --overview-chunk N   overview bins per stored row (default 256)
  --overview-sv BP     the excursion length that marks a bin variant
                       (default 50)
  --forward-only       sample only the forward orientation of each path
  --page-size BYTES    SQLite page size of index.db, a power of two from 512
                       to 65536 (default 65536: one page per 64 KiB block of
                       a remote reader, and shallower b-trees)
  --overwrite          replace index.db if it exists
  --threads N          walker threads for the GBZ route
";

const ANCHOR_RULE: &str = "HaplotypeAnchors names, per indexed reference path and multiple k of the spacing, the path's first node for k = 0 and for k >= 1 the node with the most GBWT positions among those overlapping [k * spacing - spacing / 2, k * spacing), the last on a tie; HaplotypeSamples holds every GBWT position at those nodes, in both orientations";

struct Args {
    gbz: Option<String>,
    db: Option<String>,
    output: String,
    overwrite: bool,
    interval: usize,
    anchor_spacing: usize,
    anchor_sample: Option<String>,
    reference_interval: usize,
    stray_context: usize,
    stray_bin: usize,
    stray_bound: usize,
    stray_gap: usize,
    forward_only: bool,
    page_size: usize,
    overview_bin: usize,
    overview_chunk: usize,
    overview_sv: usize,
    threads: usize,
}

impl Args {
    fn stray_options(&self) -> strays::Options {
        strays::Options {
            context: self.stray_context,
            bin: self.stray_bin,
            bound: self.stray_bound,
            gap: self.stray_gap,
            part_bytes: encode::part_limit(self.page_size) / 2,
        }
    }

    fn overview_options(&self) -> overview::Options {
        overview::Options {
            bin: self.overview_bin,
            chunk: self.overview_chunk,
            sv_bp: self.overview_sv,
            part_bytes: encode::part_limit(self.page_size),
        }
    }
}

fn parse_args() -> Args {
    let mut positional = Vec::new();
    let mut interval = 4096;
    let mut anchor_spacing = 32768;
    let mut anchor_sample = None;
    let mut reference_interval = None;
    let mut stray_context = 1000;
    let mut stray_bin = 16384;
    let mut stray_bound = 32768;
    let mut stray_gap = 1024;
    let mut forward_only = false;
    let mut page_size: usize = 65536;
    let mut overview_bin = 4096;
    let mut overview_chunk = 256;
    let mut overview_sv = 50;
    let mut from_db = false;
    let mut overwrite = false;
    let mut threads = thread::available_parallelism()
        .map(|n| n.get())
        .unwrap_or(1);
    let mut iter = env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--interval" => {
                let value = iter.next().unwrap_or_default();
                interval = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --interval: {}", value);
                    process::exit(1);
                });
            }
            "--anchor-spacing" => {
                let value = iter.next().unwrap_or_default();
                anchor_spacing = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --anchor-spacing: {}", value);
                    process::exit(1);
                });
            }
            "--anchor-sample" => {
                anchor_sample = Some(iter.next().unwrap_or_else(|| {
                    eprintln!("--anchor-sample needs a sample name");
                    process::exit(1);
                }));
            }
            "--reference-interval" => {
                let value = iter.next().unwrap_or_default();
                reference_interval = Some(value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --reference-interval: {}", value);
                    process::exit(1);
                }));
            }
            "--stray-context" | "--stray-bin" | "--stray-bound" | "--stray-gap"
            | "--overview-bin" | "--overview-chunk" | "--overview-sv" => {
                let value = iter.next().unwrap_or_default();
                let parsed = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid {}: {}", arg, value);
                    process::exit(1);
                });
                match arg.as_str() {
                    "--stray-context" => stray_context = parsed,
                    "--stray-bin" => stray_bin = parsed,
                    "--stray-bound" => stray_bound = parsed,
                    "--stray-gap" => stray_gap = parsed,
                    "--overview-bin" => overview_bin = parsed,
                    "--overview-chunk" => overview_chunk = parsed,
                    _ => overview_sv = parsed,
                }
            }
            "--threads" => {
                let value = iter.next().unwrap_or_default();
                threads = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --threads: {}", value);
                    process::exit(1);
                });
            }
            "--page-size" => {
                let value = iter.next().unwrap_or_default();
                page_size = value.parse().unwrap_or(0);
                if !(512..=65536).contains(&page_size) || !page_size.is_power_of_two() {
                    eprintln!("Invalid --page-size: {}", value);
                    process::exit(1);
                }
            }
            "--forward-only" => forward_only = true,
            "--from-db" => from_db = true,
            "--overwrite" => overwrite = true,
            "-h" | "--help" => {
                eprint!("{}", USAGE);
                process::exit(0);
            }
            _ if arg.starts_with("--") => {
                eprintln!("Unknown option {}", arg);
                eprint!("{}", USAGE);
                process::exit(1);
            }
            _ => positional.push(arg),
        }
    }
    let valid = if from_db {
        positional.len() == 2
    } else {
        positional.len() == 2 || positional.len() == 3
    };
    let reference_interval = reference_interval.unwrap_or(interval);
    if !valid
        || interval == 0
        || reference_interval == 0
        || threads == 0
        || stray_bin == 0
        || overview_chunk == 0
        || overview_sv == 0
        || overview_bin > u16::MAX as usize
    {
        eprint!("{}", USAGE);
        process::exit(1);
    }
    let output = positional.pop().unwrap();
    let (gbz, db) = if from_db {
        (None, positional.pop())
    } else {
        let gbz = positional.remove(0);
        (Some(gbz), positional.pop())
    };
    Args {
        gbz,
        db,
        output,
        overwrite,
        interval,
        anchor_spacing,
        anchor_sample,
        reference_interval,
        stray_context,
        stray_bin,
        stray_bound,
        stray_gap,
        forward_only,
        page_size,
        overview_bin,
        overview_chunk,
        overview_sv,
        threads,
    }
}

#[derive(Clone, Copy)]
struct Sample {
    node_handle: u32,
    node_offset: u32,
    path_handle: u32,
    orientation: u8,
    path_offset: u32,
}

#[derive(Clone, Copy)]
struct Anchor {
    path_handle: u32,
    anchor_offset: u32,
    node_handle: u32,
    path_offset: u32,
}

fn fit_u32<T: TryInto<u32> + Copy + std::fmt::Display>(field: &str, value: T) -> u32 {
    value.try_into().unwrap_or_else(|_| {
        eprintln!("{} {} does not fit the index's 32-bit column", field, value);
        process::exit(1);
    })
}

struct NodeSet {
    words: Vec<u64>,
}

impl NodeSet {
    fn with_capacity(max_node_id: usize) -> Self {
        NodeSet {
            words: vec![0; max_node_id / 64 + 1],
        }
    }

    fn insert(&mut self, node_id: usize) {
        self.words[node_id / 64] |= 1 << (node_id % 64);
    }

    fn contains(&self, node_id: usize) -> bool {
        self.words
            .get(node_id / 64)
            .map_or(false, |word| word & (1 << (node_id % 64)) != 0)
    }

    fn len(&self) -> usize {
        self.words
            .iter()
            .map(|word| word.count_ones() as usize)
            .sum()
    }

    #[cfg(test)]
    fn ids(&self) -> Vec<usize> {
        (0..self.words.len() * 64)
            .filter(|&id| self.contains(id))
            .collect()
    }
}

trait PathSource {
    fn start(&self, path_handle: usize, orientation: Orientation) -> Option<Pos>;
    fn step(&self, pos: Pos) -> (usize, Option<Pos>);
    fn node_len(&self, handle: usize) -> usize;
    fn visits(&self, handle: usize) -> usize;
    fn indexed_paths(&self, sample: Option<&str>) -> Vec<usize>;
    fn max_node_id(&self) -> usize;
    fn successors(&self, handle: usize) -> Vec<usize>;
    fn sample_of(&self, path_handle: usize) -> String;
    fn haplotype_of(&self, path_handle: usize) -> usize;
    fn path_count(&self) -> usize;
}

struct GbzSource<'a> {
    graph: &'a GBZ,
}

impl PathSource for GbzSource<'_> {
    fn start(&self, path_handle: usize, orientation: Orientation) -> Option<Pos> {
        let index: &GBWT = self.graph.as_ref();
        index.start(support::encode_path(path_handle, orientation))
    }

    fn step(&self, pos: Pos) -> (usize, Option<Pos>) {
        let node_len = self.graph.sequence_len(support::node_id(pos.node)).unwrap();
        let index: &GBWT = self.graph.as_ref();
        let bwt: &BWT = index.as_ref();
        let record = bwt.record(index.node_to_record(pos.node)).unwrap();
        (node_len, record.lf(pos.offset))
    }

    fn node_len(&self, handle: usize) -> usize {
        self.graph.sequence_len(support::node_id(handle)).unwrap()
    }

    fn visits(&self, handle: usize) -> usize {
        let index: &GBWT = self.graph.as_ref();
        let bwt: &BWT = index.as_ref();
        bwt.record(index.node_to_record(handle))
            .map_or(0, |record| record.len())
    }

    fn indexed_paths(&self, sample: Option<&str>) -> Vec<usize> {
        let reference_samples: BTreeSet<usize> =
            self.graph.reference_sample_ids(true).into_iter().collect();
        match self.graph.metadata() {
            Some(metadata) => metadata
                .path_iter()
                .enumerate()
                .filter(|(_, name)| reference_samples.contains(&name.sample()))
                .filter(|(_, name)| {
                    sample.map_or(true, |s| metadata.sample_name(name.sample()) == s)
                })
                .map(|(handle, _)| handle)
                .collect(),
            None => Vec::new(),
        }
    }

    fn max_node_id(&self) -> usize {
        self.graph.max_node()
    }

    fn successors(&self, handle: usize) -> Vec<usize> {
        self.graph
            .successors(support::node_id(handle), support::node_orientation(handle))
            .map_or_else(Vec::new, |edges| {
                edges.map(|(id, o)| support::encode_node(id, o)).collect()
            })
    }

    fn sample_of(&self, path_handle: usize) -> String {
        let metadata = self.graph.metadata().unwrap();
        metadata.sample_name(metadata.path(path_handle).unwrap().sample())
    }

    fn haplotype_of(&self, path_handle: usize) -> usize {
        self.graph
            .metadata()
            .and_then(|m| m.path(path_handle))
            .map_or(0, |name| name.phase())
    }

    fn path_count(&self) -> usize {
        self.graph.metadata().map_or(0, |m| m.paths())
    }
}

struct DbSource<'a> {
    interface: std::cell::RefCell<GraphInterface<'a>>,
    paths: usize,
    max_node_id: usize,
}

impl PathSource for DbSource<'_> {
    fn start(&self, path_handle: usize, orientation: Orientation) -> Option<Pos> {
        let path = self.interface.borrow_mut().get_path(path_handle).unwrap()?;
        Some(if orientation == Orientation::Forward {
            path.fw_start
        } else {
            path.rev_start
        })
    }

    fn step(&self, pos: Pos) -> (usize, Option<Pos>) {
        let record = self
            .interface
            .borrow_mut()
            .get_record(pos.node)
            .unwrap()
            .unwrap();
        (
            record.sequence_len(),
            record.to_gbwt_record().lf(pos.offset),
        )
    }

    fn node_len(&self, handle: usize) -> usize {
        self.interface
            .borrow_mut()
            .get_record(handle)
            .unwrap()
            .unwrap()
            .sequence_len()
    }

    fn visits(&self, handle: usize) -> usize {
        self.interface
            .borrow_mut()
            .get_record(handle)
            .unwrap()
            .map_or(0, |record| record.to_gbwt_record().len())
    }

    fn indexed_paths(&self, sample: Option<&str>) -> Vec<usize> {
        (0..self.paths)
            .filter(|&handle| {
                self.interface
                    .borrow_mut()
                    .get_path(handle)
                    .unwrap()
                    .map_or(false, |path| {
                        path.is_indexed && sample.map_or(true, |s| path.name.sample == s)
                    })
            })
            .collect()
    }

    fn max_node_id(&self) -> usize {
        self.max_node_id
    }

    fn successors(&self, handle: usize) -> Vec<usize> {
        self.interface
            .borrow_mut()
            .get_record(handle)
            .unwrap()
            .map_or_else(Vec::new, |record| record.successors().collect())
    }

    fn sample_of(&self, path_handle: usize) -> String {
        self.interface
            .borrow_mut()
            .get_path(path_handle)
            .unwrap()
            .map(|path| path.name.sample)
            .unwrap_or_default()
    }

    fn haplotype_of(&self, path_handle: usize) -> usize {
        self.interface
            .borrow_mut()
            .get_path(path_handle)
            .unwrap()
            .map_or(0, |path| path.name.haplotype)
    }

    fn path_count(&self) -> usize {
        self.paths
    }
}

// One anchor per multiple k of the spacing along a reference path: the first
// node for k = 0, and for k >= 1 the node with the most GBWT positions among
// those overlapping [k * spacing - spacing / 2, k * spacing), the last on a
// tie. A window starting at or past k * spacing has that anchor before it, and
// most haplotypes of the region visit it, where the node that happens to
// contain the multiple can be a rare allele. Most reference nodes tie, so the
// last one keeps the anchor close to the multiple and the walk from it short.
fn mark_anchors(
    source: &dyn PathSource,
    path_handle: usize,
    spacing: usize,
    anchors: &mut Vec<Anchor>,
) {
    let half = spacing / 2;
    let mut pos = source.start(path_handle, Orientation::Forward);
    let mut offset = 0;
    let mut k = 1;
    let mut best: Option<(usize, Pos, usize)> = None;
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (node_len, next) = source.step(current);
        if offset == 0 {
            anchors.push(Anchor {
                path_handle: path_handle as u32,
                anchor_offset: 0,
                node_handle: fit_u32("node_handle", current.node),
                path_offset: 0,
            });
        }
        let end = offset + node_len;
        let mut visits = None;
        loop {
            let target = k * spacing;
            if end <= target - half {
                break;
            }
            let count = *visits.get_or_insert_with(|| source.visits(current.node));
            if best.map_or(true, |(most, _, _)| count >= most) {
                best = Some((count, current, offset));
            }
            if end >= target {
                if let Some((_, chosen, chosen_offset)) = best.take() {
                    anchors.push(Anchor {
                        path_handle: path_handle as u32,
                        anchor_offset: fit_u32("anchor_offset", target),
                        node_handle: fit_u32("node_handle", chosen.node),
                        path_offset: fit_u32("path_offset", chosen_offset),
                    });
                }
                k += 1;
            } else {
                break;
            }
        }
        offset = end;
        pos = next;
    }
}

fn anchor_set(source: &dyn PathSource, anchors: &[Anchor]) -> NodeSet {
    let mut nodes = NodeSet::with_capacity(source.max_node_id());
    for anchor in anchors {
        nodes.insert(support::node_id(anchor.node_handle as usize));
    }
    nodes
}

fn anchor_rows(source: &dyn PathSource, reference_paths: &[usize], spacing: usize) -> Vec<Anchor> {
    let mut anchors = Vec::new();
    if spacing > 0 {
        for &path_handle in reference_paths {
            mark_anchors(source, path_handle, spacing, &mut anchors);
        }
    }
    anchors
}

fn walk(
    source: &dyn PathSource,
    path_handle: usize,
    orientation: Orientation,
    interval: usize,
    anchors: &NodeSet,
    samples: &mut Vec<Sample>,
) -> usize {
    let mut pos = source.start(path_handle, orientation);
    let mut offset = 0;
    let mut next_sample = 0;
    let mut last: Option<(Pos, usize)> = None;
    let first = samples.len();
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (node_len, next) = source.step(current);
        if offset >= next_sample || anchors.contains(support::node_id(current.node)) {
            samples.push(Sample {
                node_handle: fit_u32("node_handle", current.node),
                node_offset: fit_u32("node_offset", current.offset),
                path_handle: path_handle as u32,
                orientation: orientation as u8,
                path_offset: fit_u32("path_offset", offset),
            });
            next_sample = offset + interval;
            last = None;
        } else {
            last = Some((current, offset));
        }
        offset += node_len;
        pos = next;
    }
    if let Some((end, end_offset)) = last {
        samples.push(Sample {
            node_handle: fit_u32("node_handle", end.node),
            node_offset: fit_u32("node_offset", end.offset),
            path_handle: path_handle as u32,
            orientation: orientation as u8,
            path_offset: fit_u32("path_offset", end_offset),
        });
    }
    let length = offset;
    if orientation == Orientation::Reverse {
        for sample in samples[first..].iter_mut() {
            let node_len = source.node_len(sample.node_handle as usize);
            sample.path_offset = fit_u32(
                "path_offset",
                length - sample.path_offset as usize - node_len,
            );
        }
    }
    length
}

fn orientations(args: &Args) -> Vec<Orientation> {
    if args.forward_only {
        vec![Orientation::Forward]
    } else {
        vec![Orientation::Forward, Orientation::Reverse]
    }
}

fn walk_paths(
    source: &dyn PathSource,
    handles: std::ops::Range<usize>,
    args: &Args,
    anchors: &Anchors,
    label: &str,
) -> (Vec<Sample>, Vec<(usize, usize)>) {
    let mut samples = Vec::new();
    let mut lengths = Vec::new();
    let started = Instant::now();
    let total = handles.len();
    let mut walked_bp: usize = 0;
    for (done, path_handle) in handles.enumerate() {
        let mut length = 0;
        let interval = if anchors.reference_paths.contains(&path_handle) {
            args.reference_interval
        } else {
            args.interval
        };
        for &orientation in orientations(args).iter() {
            length = walk(
                source,
                path_handle,
                orientation,
                interval,
                &anchors.nodes,
                &mut samples,
            );
        }
        walked_bp += length;
        lengths.push((path_handle, length));
        if (done + 1) % 500 == 0 || done + 1 == total {
            eprintln!(
                "{}: {} / {} paths, {:.2} Gbp, {} samples, {:.0} s",
                label,
                done + 1,
                total,
                walked_bp as f64 / 1e9,
                samples.len(),
                started.elapsed().as_secs_f64()
            );
        }
    }
    (samples, lengths)
}

struct Anchors {
    nodes: NodeSet,
    rows: Vec<Anchor>,
    reference_paths: BTreeSet<usize>,
}

fn anchored_paths(source: &dyn PathSource, args: &Args) -> Vec<usize> {
    let paths = source.indexed_paths(args.anchor_sample.as_deref());
    if let (Some(sample), true, true) = (
        &args.anchor_sample,
        paths.is_empty(),
        args.anchor_spacing > 0,
    ) {
        eprintln!("No reference path belongs to sample {}; drop --anchor-sample or name a reference sample", sample);
        process::exit(1);
    }
    paths
}

fn anchors_gbz(graph: &GBZ, args: &Args) -> Anchors {
    let started = Instant::now();
    let source = GbzSource { graph };
    let reference_paths = anchored_paths(&source, args);
    let chunk = (reference_paths.len() + args.threads - 1) / args.threads.max(1);
    let mut rows = Vec::new();
    if args.anchor_spacing > 0 && chunk > 0 {
        let marked: Vec<Vec<Anchor>> = thread::scope(|scope| {
            let workers: Vec<_> = reference_paths
                .chunks(chunk)
                .map(|handles| {
                    scope.spawn(move || {
                        let source = GbzSource { graph };
                        let mut anchors = Vec::new();
                        for &path_handle in handles {
                            mark_anchors(&source, path_handle, args.anchor_spacing, &mut anchors);
                        }
                        anchors
                    })
                })
                .collect();
            workers.into_iter().map(|w| w.join().unwrap()).collect()
        });
        rows = marked.into_iter().flatten().collect();
    }
    let nodes = anchor_set(&source, &rows);
    eprintln!(
        "Chose {} anchors on {} distinct nodes over {} reference paths at {} bp spacing in {:.0} s",
        rows.len(),
        nodes.len(),
        reference_paths.len(),
        args.anchor_spacing,
        started.elapsed().as_secs_f64()
    );
    Anchors {
        nodes,
        rows,
        reference_paths: reference_paths.into_iter().collect(),
    }
}

fn walk_gbz(
    graph: &GBZ,
    paths: usize,
    args: &Args,
    anchors: &Anchors,
) -> (Vec<Sample>, Vec<(usize, usize)>) {
    let chunk = (paths + args.threads - 1) / args.threads;
    let started = Instant::now();
    let results: Vec<(Vec<Sample>, Vec<(usize, usize)>)> = thread::scope(|scope| {
        let workers: Vec<_> = (0..args.threads)
            .map(|t| {
                let range = (t * chunk).min(paths)..((t + 1) * chunk).min(paths);
                scope.spawn(move || {
                    let source = GbzSource { graph };
                    let label = format!("thread {} (paths {}..{})", t, range.start, range.end);
                    let result = walk_paths(&source, range.clone(), args, anchors, &label);
                    eprintln!("{} done in {:.0} s", label, started.elapsed().as_secs_f64());
                    result
                })
            })
            .collect();
        workers.into_iter().map(|w| w.join().unwrap()).collect()
    });
    let mut samples = Vec::new();
    let mut lengths = Vec::new();
    for (s, l) in results {
        samples.extend(s);
        lengths.extend(l);
    }
    lengths.sort_unstable();
    (samples, lengths)
}

const SCHEMA: &str = "CREATE TABLE HaplotypeSamples (
    node_handle INTEGER NOT NULL,
    node_offset INTEGER NOT NULL,
    path_handle INTEGER NOT NULL,
    orientation INTEGER NOT NULL,
    path_offset INTEGER NOT NULL,
    PRIMARY KEY (node_handle, node_offset)
) STRICT, WITHOUT ROWID;
CREATE TABLE HaplotypeLengths (
    path_handle INTEGER PRIMARY KEY,
    length INTEGER NOT NULL
) STRICT;
CREATE TABLE Tags (
    key TEXT PRIMARY KEY,
    value TEXT NOT NULL
) STRICT;
CREATE TABLE HaplotypeAnchors (
    path_handle INTEGER NOT NULL,
    anchor_offset INTEGER NOT NULL,
    part INTEGER NOT NULL,
    node_handle INTEGER NOT NULL,
    path_offset INTEGER NOT NULL,
    visits BLOB NOT NULL,
    PRIMARY KEY (path_handle, anchor_offset, part)
) STRICT;
CREATE TABLE HaplotypeBins (
    reference_handle INTEGER NOT NULL,
    bin INTEGER NOT NULL,
    part INTEGER NOT NULL,
    nodes BLOB NOT NULL,
    strays BLOB NOT NULL,
    PRIMARY KEY (reference_handle, bin, part)
) STRICT;
CREATE TABLE HaplotypeOverviewRows (
    row INTEGER PRIMARY KEY,
    sample TEXT NOT NULL,
    haplotype INTEGER NOT NULL
) STRICT;
CREATE TABLE HaplotypeOverviewBins (
    reference_handle INTEGER NOT NULL,
    level INTEGER NOT NULL,
    chunk INTEGER NOT NULL,
    part INTEGER NOT NULL,
    bins BLOB NOT NULL,
    PRIMARY KEY (reference_handle, level, chunk, part)
) STRICT;
CREATE TABLE HaplotypeOverviewClasses (
    reference_handle INTEGER NOT NULL,
    level INTEGER NOT NULL,
    chunk INTEGER NOT NULL,
    part INTEGER NOT NULL,
    classes BLOB NOT NULL,
    PRIMARY KEY (reference_handle, level, chunk, part)
) STRICT;";

// The layout of the tables, which the reader tests for. Format 3 stores the
// samples in their key b-tree, each anchor with the visits through its node
// and each bin with its stray rows, so a window's rows come in a few reads.
// The tables with blobs keep SQLite's rowid and a key index, since an index
// b-tree copies whole rows into its interior pages and a blob would leave
// each page a few cells; their rows go in key order, so a key range is one
// run of the table.
const FORMAT: &str = "3";

// The forward samples at both handles of an anchor's node, from the samples
// sorted by (node_handle, node_offset).
fn visits_at(samples: &[Sample], node_handle: u32) -> Vec<Sample> {
    let id = support::node_id(node_handle as usize) as u32;
    let mut visits = Vec::new();
    for handle in [2 * id, 2 * id + 1] {
        let from = samples.partition_point(|s| s.node_handle < handle);
        visits.extend(
            samples[from..]
                .iter()
                .take_while(|s| s.node_handle == handle)
                .filter(|s| s.orientation == Orientation::Forward as u8)
                .copied(),
        );
    }
    visits
}

fn write(
    target: &str,
    mut samples: Vec<Sample>,
    lengths: &[(usize, usize)],
    paths: usize,
    nodes: usize,
    anchors: &Anchors,
    strays: &strays::Output,
    overview: &overview::Output,
    args: &Args,
) {
    let started = Instant::now();
    samples.sort_unstable_by_key(|s| (s.node_handle, s.node_offset));
    eprintln!(
        "Sorted {} samples in {:.0} s",
        samples.len(),
        started.elapsed().as_secs_f64()
    );
    let partial = format!("{}.partial", target);
    let _ = std::fs::remove_file(&partial);
    let mut connection = Connection::open(&partial).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", partial, e);
        process::exit(1);
    });
    connection
        .execute_batch(&format!(
            "PRAGMA page_size = {}; PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF;",
            args.page_size
        ))
        .unwrap();
    connection.execute_batch(SCHEMA).unwrap();
    let limit = encode::part_limit(args.page_size);
    let transaction = connection.transaction().unwrap();
    let mut anchor_parts = 0;
    let mut bin_parts = 0;
    let mut bin_bytes = 0;
    {
        let mut write_sample = transaction
            .prepare("INSERT INTO HaplotypeSamples(node_handle, node_offset, path_handle, orientation, path_offset) VALUES (?1, ?2, ?3, ?4, ?5)")
            .unwrap();
        for s in samples.iter() {
            write_sample
                .execute(params![
                    s.node_handle as i64,
                    s.node_offset as i64,
                    s.path_handle as i64,
                    s.orientation as i64,
                    s.path_offset as i64
                ])
                .unwrap();
        }
        let mut write_length = transaction
            .prepare("INSERT INTO HaplotypeLengths(path_handle, length) VALUES (?1, ?2)")
            .unwrap();
        for &(handle, length) in lengths {
            write_length
                .execute(params![handle as i64, length as i64])
                .unwrap();
        }
        let mut write_anchor = transaction
            .prepare("INSERT INTO HaplotypeAnchors(path_handle, anchor_offset, part, node_handle, path_offset, visits) VALUES (?1, ?2, ?3, ?4, ?5, ?6)")
            .unwrap();
        let mut rows: Vec<&Anchor> = anchors.rows.iter().collect();
        rows.sort_unstable_by_key(|a| (a.path_handle, a.anchor_offset));
        for a in rows {
            let mut parts = encode::encode_visits(&visits_at(&samples, a.node_handle), limit);
            if parts.is_empty() {
                parts.push(Vec::new());
            }
            for (part, visits) in parts.iter().enumerate() {
                write_anchor
                    .execute(params![
                        a.path_handle as i64,
                        a.anchor_offset as i64,
                        part as i64,
                        a.node_handle as i64,
                        a.path_offset as i64,
                        visits
                    ])
                    .unwrap();
                anchor_parts += 1;
            }
        }
        let mut write_bin = transaction.prepare("INSERT INTO HaplotypeBins(reference_handle, bin, part, nodes, strays) VALUES (?1, ?2, ?3, ?4, ?5)").unwrap();
        let mut by_bin: BTreeMap<(u32, u32), (Vec<&Vec<u8>>, Vec<strays::Stray>)> = BTreeMap::new();
        for b in strays.bins.iter() {
            by_bin
                .entry((b.reference_handle, b.bin))
                .or_default()
                .0
                .push(&b.nodes);
        }
        for r in strays.rows.iter() {
            by_bin
                .entry((r.reference_handle, r.bin))
                .or_default()
                .1
                .push(*r);
        }
        let empty = Vec::new();
        for ((reference_handle, bin), (node_parts, rows)) in by_bin {
            let stray_parts = encode::encode_strays(&rows, limit / 2);
            for part in 0..node_parts.len().max(stray_parts.len()) {
                let nodes = node_parts.get(part).copied().unwrap_or(&empty);
                let strays = stray_parts.get(part).unwrap_or(&empty);
                write_bin
                    .execute(params![
                        reference_handle as i64,
                        bin as i64,
                        part as i64,
                        nodes,
                        strays
                    ])
                    .unwrap();
                bin_parts += 1;
                bin_bytes += nodes.len() + strays.len();
            }
        }
        let mut write_row = transaction
            .prepare(
                "INSERT INTO HaplotypeOverviewRows(row, sample, haplotype) VALUES (?1, ?2, ?3)",
            )
            .unwrap();
        for (row, name) in overview.rows.iter().enumerate() {
            write_row
                .execute(params![row as i64, name.sample, name.haplotype as i64])
                .unwrap();
        }
        for (table, column, parts) in [
            ("HaplotypeOverviewBins", "bins", &overview.bins),
            ("HaplotypeOverviewClasses", "classes", &overview.classes),
        ] {
            let mut write_part = transaction
                .prepare(&format!(
                    "INSERT INTO {}(reference_handle, level, chunk, part, {}) VALUES (?1, ?2, ?3, ?4, ?5)",
                    table, column
                ))
                .unwrap();
            let mut ordered: Vec<&overview::Part> = parts.iter().collect();
            ordered.sort_unstable_by_key(|p| (p.reference_handle, p.level, p.chunk, p.part));
            for p in ordered {
                write_part
                    .execute(params![
                        p.reference_handle as i64,
                        p.level as i64,
                        p.chunk as i64,
                        p.part as i64,
                        p.bytes
                    ])
                    .unwrap();
            }
        }
        let mut write_tag = transaction
            .prepare("INSERT INTO Tags(key, value) VALUES (?1, ?2)")
            .unwrap();
        let tags = [
            ("haplotype_index_format", FORMAT.to_string()),
            (
                "haplotype_index_tool_version",
                env!("CARGO_PKG_VERSION").to_string(),
            ),
            ("haplotype_index_interval", args.interval.to_string()),
            (
                "haplotype_index_reference_interval",
                args.reference_interval.to_string(),
            ),
            (
                "haplotype_index_orientations",
                if args.forward_only { "forward" } else { "both" }.to_string(),
            ),
            ("haplotype_index_paths", paths.to_string()),
            ("haplotype_index_nodes", nodes.to_string()),
        ];
        for (key, value) in tags.iter() {
            write_tag.execute(params![key, value]).unwrap();
        }
        if args.anchor_spacing > 0 {
            let anchor_tags = [
                (
                    "haplotype_index_anchor_spacing",
                    args.anchor_spacing.to_string(),
                ),
                ("haplotype_index_anchor_rule", ANCHOR_RULE.to_string()),
                (
                    "haplotype_index_anchor_paths",
                    anchors.reference_paths.len().to_string(),
                ),
                (
                    "haplotype_index_anchor_nodes",
                    anchors.nodes.len().to_string(),
                ),
                (
                    "haplotype_index_anchor_sample",
                    args.anchor_sample.clone().unwrap_or_default(),
                ),
            ];
            for (key, value) in anchor_tags.iter() {
                write_tag.execute(params![key, value]).unwrap();
            }
        }
        if !strays.samples.is_empty() {
            let stray_tags = [
                ("haplotype_index_stray_format", "2".to_string()),
                (
                    "haplotype_index_stray_context",
                    args.stray_context.to_string(),
                ),
                ("haplotype_index_stray_bin", args.stray_bin.to_string()),
                ("haplotype_index_stray_bound", args.stray_bound.to_string()),
                ("haplotype_index_stray_gap", args.stray_gap.to_string()),
                ("haplotype_index_stray_samples", strays.samples.join(",")),
                ("haplotype_index_stray_rows", strays.rows.len().to_string()),
                ("haplotype_index_stray_bin_parts", bin_parts.to_string()),
                (
                    "haplotype_index_stray_snarls",
                    if strays.snarls_modeled {
                        "modeled"
                    } else {
                        "none"
                    }
                    .to_string(),
                ),
                (
                    "haplotype_index_stray_snarl_nodes",
                    strays::SNARL_NODES.to_string(),
                ),
                (
                    "haplotype_index_stray_chain_links",
                    strays.chain_links.to_string(),
                ),
                ("haplotype_index_stray_rule", strays::RULE.to_string()),
            ];
            for (key, value) in stray_tags.iter() {
                write_tag.execute(params![key, value]).unwrap();
            }
        }
        if !overview.samples.is_empty() {
            let overview_tags = [
                (
                    "haplotype_index_overview_format",
                    overview::FORMAT.to_string(),
                ),
                (
                    "haplotype_index_overview_bin",
                    args.overview_bin.to_string(),
                ),
                (
                    "haplotype_index_overview_chunk",
                    args.overview_chunk.to_string(),
                ),
                ("haplotype_index_overview_sv", args.overview_sv.to_string()),
                (
                    "haplotype_index_overview_levels",
                    overview.levels.to_string(),
                ),
                (
                    "haplotype_index_overview_rows",
                    overview.rows.len().to_string(),
                ),
                (
                    "haplotype_index_overview_samples",
                    overview.samples.join(","),
                ),
            ];
            for (key, value) in overview_tags.iter() {
                write_tag.execute(params![key, value]).unwrap();
            }
        }
    }
    transaction.commit().unwrap();
    if let Err((_, e)) = connection.close() {
        eprintln!("Cannot close {}: {}", partial, e);
        process::exit(1);
    }
    std::fs::rename(&partial, target).unwrap_or_else(|e| {
        eprintln!("Cannot move {} to {}: {}", partial, target, e);
        process::exit(1);
    });
    eprintln!(
        "Wrote {} samples, {} anchors in {} parts, {} stray rows and {} bin parts ({} bytes), an overview of {} haplotypes in {} parts ({} bytes) for {} paths to {} in {:.0} s",
        samples.len(),
        anchors.rows.len(),
        anchor_parts,
        strays.rows.len(),
        bin_parts,
        bin_bytes,
        overview.rows.len(),
        overview.bins.len() + overview.classes.len(),
        overview
            .bins
            .iter()
            .chain(overview.classes.iter())
            .map(|p| p.bytes.len())
            .sum::<usize>(),
        paths,
        target,
        started.elapsed().as_secs_f64()
    );
}

fn open_read_only(db: &str) -> Connection {
    Connection::open_with_flags(db, OpenFlags::SQLITE_OPEN_READ_ONLY).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", db, e);
        process::exit(1);
    })
}

fn count_from_db(db: &str, key: &str) -> usize {
    let connection = open_read_only(db);
    let value: String = connection
        .query_row(
            "SELECT value FROM Tags WHERE key = ?1",
            params![key],
            |row| row.get(0),
        )
        .unwrap_or_else(|_| "0".to_string());
    value.parse().unwrap_or(0)
}

fn max_node_id_from_db(db: &str) -> usize {
    let connection = open_read_only(db);
    let handle: i64 = connection
        .query_row("SELECT max(handle) FROM Nodes", [], |row| row.get(0))
        .unwrap_or(0);
    support::node_id(handle as usize)
}

// The links between the boundary nodes of top-level snarls, as (handle, next).
fn chain_links(db: &str) -> Vec<(usize, usize)> {
    let connection = open_read_only(db);
    let mut statement = connection
        .prepare("SELECT handle, next FROM Nodes WHERE next IS NOT NULL")
        .unwrap();
    let links = statement
        .query_map([], |row| {
            Ok((
                row.get::<_, i64>(0)? as usize,
                row.get::<_, i64>(1)? as usize,
            ))
        })
        .unwrap();
    links.map(|link| link.unwrap()).collect()
}

fn walk_db(
    db: &str,
    args: &Args,
) -> (
    Vec<Sample>,
    Vec<(usize, usize)>,
    Anchors,
    strays::Output,
    overview::Output,
) {
    let paths = count_from_db(db, "paths");
    let max_node_id = max_node_id_from_db(db);
    let database = GBZBase::open(db).unwrap_or_else(|e| {
        eprintln!("Cannot open {} as a GBZ-base: {}", db, e);
        process::exit(1);
    });
    let interface = GraphInterface::new(&database).unwrap();
    let source = DbSource {
        interface: std::cell::RefCell::new(interface),
        paths,
        max_node_id,
    };
    let reference_paths = anchored_paths(&source, args);
    let rows = anchor_rows(&source, &reference_paths, args.anchor_spacing);
    let nodes = anchor_set(&source, &rows);
    eprintln!(
        "Chose {} anchors on {} distinct nodes over {} reference paths at {} bp spacing",
        rows.len(),
        nodes.len(),
        reference_paths.len(),
        args.anchor_spacing
    );
    let anchors = Anchors {
        nodes,
        rows,
        reference_paths: reference_paths.into_iter().collect(),
    };
    let (samples, lengths) = walk_paths(&source, 0..paths, args, &anchors, "database walk");
    let strays = if args.stray_context > 0 && args.anchor_spacing > 0 {
        let snarls =
            strays::Snarls::find(&source, &chain_links(db), count_from_db(db, "chain_links"));
        strays::strays(
            &strays::Serial(&source),
            &anchors,
            args.anchor_spacing,
            args.stray_options(),
            &snarls,
        )
    } else {
        strays::Output::empty()
    };
    let overview = if args.anchor_spacing > 0 {
        overview::overview(
            &strays::Serial(&source),
            &anchors,
            &lengths,
            args.overview_options(),
        )
    } else {
        overview::Output::empty()
    };
    (samples, lengths, anchors, strays, overview)
}

fn main() {
    let args = parse_args();
    if std::path::Path::new(&args.output).exists() && !args.overwrite {
        eprintln!("{} exists; pass --overwrite to replace it", args.output);
        process::exit(1);
    }
    let (samples, lengths, paths, nodes, anchors, strays, overview) = match &args.gbz {
        Some(gbz) => {
            let started = Instant::now();
            let graph: GBZ = serialize::load_from(gbz).unwrap_or_else(|e| {
                eprintln!("Cannot load {}: {}", gbz, e);
                process::exit(1);
            });
            let paths = graph.metadata().map(|m| m.paths()).unwrap_or(0);
            if paths == 0 {
                eprintln!("The GBZ has no path metadata");
                process::exit(1);
            }
            eprintln!(
                "Loaded {} with {} paths in {:.0} s",
                gbz,
                paths,
                started.elapsed().as_secs_f64()
            );
            let nodes = graph.nodes();
            if let Some(db) = &args.db {
                let db_paths = count_from_db(db, "paths");
                if db_paths != paths {
                    eprintln!("{} has {} paths but {} has {}", gbz, paths, db, db_paths);
                    process::exit(1);
                }
                let db_nodes = count_from_db(db, "nodes");
                if db_nodes != nodes {
                    eprintln!("{} has {} nodes but {} has {}", gbz, nodes, db, db_nodes);
                    process::exit(1);
                }
            }
            let anchors = anchors_gbz(&graph, &args);
            let (samples, lengths) = walk_gbz(&graph, paths, &args, &anchors);
            let source = GbzSource { graph: &graph };
            let strays = if args.stray_context > 0 && args.anchor_spacing > 0 {
                let snarls = match &args.db {
                    Some(db) => strays::Snarls::find(
                        &source,
                        &chain_links(db),
                        count_from_db(db, "chain_links"),
                    ),
                    None => {
                        eprintln!("No graph.gbz.db given: the stray rows leave snarls out, and a query that fills one identifies every walk");
                        strays::Snarls::none()
                    }
                };
                strays::strays(
                    &strays::Threads {
                        source: &source,
                        threads: args.threads,
                    },
                    &anchors,
                    args.anchor_spacing,
                    args.stray_options(),
                    &snarls,
                )
            } else {
                strays::Output::empty()
            };
            let overview = if args.anchor_spacing > 0 {
                overview::overview(
                    &strays::Threads {
                        source: &source,
                        threads: args.threads,
                    },
                    &anchors,
                    &lengths,
                    args.overview_options(),
                )
            } else {
                overview::Output::empty()
            };
            (samples, lengths, paths, nodes, anchors, strays, overview)
        }
        None => {
            let db = args.db.as_ref().unwrap();
            let paths = count_from_db(db, "paths");
            let nodes = count_from_db(db, "nodes");
            let (samples, lengths, anchors, strays, overview) = walk_db(db, &args);
            (samples, lengths, paths, nodes, anchors, strays, overview)
        }
    };
    write(
        &args.output,
        samples,
        &lengths,
        paths,
        nodes,
        &anchors,
        &strays,
        &overview,
        &args,
    );
}

#[cfg(test)]
mod tests {
    use super::*;

    fn split_contig() -> String {
        format!(
            "{}/test/data/split-contig.gbz.db",
            env!("CARGO_MANIFEST_DIR")
        )
    }

    fn args(interval: usize, anchor_spacing: usize) -> Args {
        Args {
            gbz: None,
            db: Some(split_contig()),
            output: String::new(),
            overwrite: false,
            interval,
            anchor_spacing,
            anchor_sample: None,
            reference_interval: interval,
            stray_context: 1000,
            stray_bin: 16384,
            stray_bound: 32768,
            stray_gap: 1024,
            forward_only: false,
            page_size: 4096,
            overview_bin: 100,
            overview_chunk: 4,
            overview_sv: 50,
            threads: 1,
        }
    }

    fn forward_nodes(source: &dyn PathSource, path_handle: usize) -> Vec<(usize, usize, usize)> {
        let mut nodes = Vec::new();
        let mut pos = source.start(path_handle, Orientation::Forward);
        let mut offset = 0;
        while let Some(current) = pos {
            if current.node == ENDMARKER {
                break;
            }
            let (len, next) = source.step(current);
            nodes.push((current.node, offset, len));
            offset += len;
            pos = next;
        }
        nodes
    }

    fn with_source<T>(f: impl FnOnce(&DbSource) -> T) -> T {
        let db = split_contig();
        let paths = count_from_db(&db, "paths");
        let max_node_id = max_node_id_from_db(&db);
        let database = GBZBase::open(&db).unwrap();
        let interface = GraphInterface::new(&database).unwrap();
        let source = DbSource {
            interface: std::cell::RefCell::new(interface),
            paths,
            max_node_id,
        };
        f(&source)
    }

    #[test]
    fn every_reference_path_gets_its_first_node_and_the_most_visited_node_nearest_each_multiple() {
        let spacing = 300;
        let (_, _, anchors, _, _) = walk_db(&split_contig(), &args(200, spacing));
        assert_eq!(anchors.reference_paths.len(), 2);
        let checked = with_source(|source| {
            let mut checked = 0;
            for &path_handle in source.indexed_paths(None).iter() {
                let nodes = forward_nodes(source, path_handle);
                let length: usize = nodes.iter().map(|n| n.2).sum();
                let first = anchors
                    .rows
                    .iter()
                    .find(|a| a.path_handle as usize == path_handle && a.anchor_offset == 0)
                    .unwrap();
                assert_eq!(
                    (first.node_handle as usize, first.path_offset),
                    (nodes[0].0, 0)
                );
                let mut k = 1;
                while k * spacing <= length {
                    let target = k * spacing;
                    let overlapping: Vec<&(usize, usize, usize)> = nodes
                        .iter()
                        .filter(|(_, start, len)| {
                            *start < target && start + len > target - spacing / 2
                        })
                        .collect();
                    let most = overlapping
                        .iter()
                        .map(|(handle, _, _)| source.visits(*handle))
                        .max()
                        .unwrap();
                    let chosen = anchors
                        .rows
                        .iter()
                        .find(|a| {
                            a.path_handle as usize == path_handle
                                && a.anchor_offset as usize == target
                        })
                        .unwrap();
                    let (handle, start, _) = overlapping
                        .iter()
                        .find(|(handle, _, _)| *handle == chosen.node_handle as usize)
                        .unwrap();
                    assert_eq!(chosen.path_offset as usize, *start);
                    assert_eq!(source.visits(*handle), most);
                    assert!(overlapping
                        .iter()
                        .rev()
                        .take_while(|(h, _, _)| h != handle)
                        .all(|(h, _, _)| source.visits(*h) < most));
                    checked += 1;
                    k += 1;
                }
            }
            checked
        });
        assert!(checked >= 2);
        let distinct: BTreeSet<usize> = anchors
            .rows
            .iter()
            .map(|a| support::node_id(a.node_handle as usize))
            .collect();
        assert_eq!(anchors.nodes.len(), distinct.len());
    }

    #[test]
    fn every_visit_through_an_anchor_node_is_sampled_in_both_orientations() {
        let (samples, lengths, anchors, _, _) = walk_db(&split_contig(), &args(200, 300));
        assert_eq!(lengths.len(), 6);
        let mut rows_by_handle = std::collections::BTreeMap::new();
        for sample in samples.iter() {
            if anchors
                .nodes
                .contains(support::node_id(sample.node_handle as usize))
            {
                rows_by_handle
                    .entry(sample.node_handle as usize)
                    .or_insert_with(Vec::new)
                    .push(sample.node_offset as usize);
            }
        }
        let expected_handles: BTreeSet<usize> = anchors
            .nodes
            .ids()
            .iter()
            .flat_map(|&id| [2 * id, 2 * id + 1])
            .collect();
        assert_eq!(
            rows_by_handle.keys().copied().collect::<BTreeSet<_>>(),
            expected_handles
        );
        with_source(|source| {
            for (handle, offsets) in rows_by_handle.iter_mut() {
                offsets.sort_unstable();
                assert_eq!(*offsets, (0..source.visits(*handle)).collect::<Vec<_>>());
            }
        });
        let mut keys: Vec<(u32, u32)> = samples
            .iter()
            .map(|s| (s.node_handle, s.node_offset))
            .collect();
        let before = keys.len();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), before);
    }

    #[test]
    fn an_anchor_sample_keeps_the_anchors_of_that_samples_paths() {
        let mut only = args(200, 300);
        only.anchor_sample = Some("GRCh38".to_string());
        let (_, _, anchors, _, _) = walk_db(&split_contig(), &only);
        let (_, _, all, _, _) = walk_db(&split_contig(), &args(200, 300));
        assert_eq!(anchors.reference_paths.len(), 2);
        assert_eq!(anchors.rows.len(), all.rows.len());
    }

    #[test]
    fn a_reference_interval_samples_only_the_reference_paths_more_densely() {
        let mut dense = args(200, 0);
        dense.reference_interval = 20;
        let (with, _, anchors, _, _) = walk_db(&split_contig(), &dense);
        let (without, _, _, _, _) = walk_db(&split_contig(), &args(200, 0));
        let per_path = |samples: &[Sample], handle: usize| -> Vec<u32> {
            let mut offsets: Vec<u32> = samples
                .iter()
                .filter(|s| {
                    s.path_handle as usize == handle && s.orientation == Orientation::Forward as u8
                })
                .map(|s| s.path_offset)
                .collect();
            offsets.sort_unstable();
            offsets
        };
        assert_eq!(anchors.reference_paths.len(), 2);
        with_source(|source| {
            for handle in 0..6 {
                if !anchors.reference_paths.contains(&handle) {
                    assert_eq!(per_path(&with, handle), per_path(&without, handle));
                    continue;
                }
                let nodes = forward_nodes(source, handle);
                let longest = nodes.iter().map(|n| n.2).max().unwrap() as u32;
                let offsets = per_path(&with, handle);
                assert!(offsets.len() > per_path(&without, handle).len());
                assert!(offsets.windows(2).all(|w| w[1] - w[0] < 20 + longest));
            }
        });
    }

    #[test]
    fn the_overview_classes_every_haplotype_in_every_bin_at_every_level() {
        let (_, _, _, _, overview) = walk_db(&split_contig(), &args(200, 300));
        let rows = overview.rows.len();
        assert!(rows >= 3);
        assert!(overview.rows.windows(2).all(|w| w[0] < w[1]));
        assert!(overview.levels >= 2);
        let bytes_per_bin = rows.div_ceil(2);
        let mut decoded = std::collections::BTreeMap::new();
        for part in &overview.classes {
            decoded
                .entry((part.reference_handle, part.level, part.chunk))
                .or_insert_with(Vec::<u8>::new)
                .extend(&part.bytes);
        }
        let mut cells = 0;
        for ((reference, level, chunk), bytes) in &decoded {
            assert_eq!(bytes.len() % bytes_per_bin, 0);
            let bins = bytes.len() / bytes_per_bin;
            assert!(
                bins <= 4
                    && (bins == 4
                        || *chunk as usize
                            == decoded
                                .keys()
                                .filter(|k| k.0 == *reference && k.1 == *level)
                                .count()
                                - 1)
            );
            for b in 0..bins {
                for row in 0..rows {
                    let byte = bytes[b * bytes_per_bin + row / 2];
                    let value = (byte >> (4 * (row % 2))) & 0xf;
                    let class = value & 3;
                    assert!(class <= overview::VARIANT);
                    assert!(class == overview::VARIANT || value >> 2 == 0);
                    let name = &overview.rows[row];
                    let own = with_source(|source| {
                        source.sample_of(*reference as usize) == name.sample
                            && source.haplotype_of(*reference as usize) == name.haplotype
                    });
                    if own {
                        assert_eq!(class, overview::REFERENCE);
                    }
                    cells += 1;
                }
            }
        }
        assert!(cells > 0);
        let mut bins_decoded = std::collections::BTreeMap::new();
        for part in &overview.bins {
            bins_decoded
                .entry((part.reference_handle, part.level, part.chunk))
                .or_insert_with(Vec::<u8>::new)
                .extend(&part.bytes);
        }
        assert_eq!(
            bins_decoded.keys().collect::<Vec<_>>(),
            decoded.keys().collect::<Vec<_>>()
        );
        for (key, bytes) in &bins_decoded {
            let values = crate::encode::tests::varints(std::slice::from_ref(bytes));
            assert_eq!(values.len() % 7, 0);
            assert_eq!(values.len() / 7, decoded[key].len() / bytes_per_bin);
            for summary in values.chunks(7) {
                assert_eq!(summary[..4].iter().sum::<usize>(), rows);
            }
        }
    }

    #[test]
    fn a_zero_spacing_marks_nothing_and_leaves_the_per_path_samples_alone() {
        let (with, _, anchors, _, _) = walk_db(&split_contig(), &args(200, 0));
        assert_eq!(anchors.nodes.len(), 0);
        assert!(anchors.rows.is_empty());
        assert_eq!(with.len(), 42);
    }
}
