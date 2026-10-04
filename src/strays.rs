// Tables HaplotypeBinNodes and HaplotypeStrays, for each reference sample with
// anchors. Together they let the keep route prove that its walks pass every
// visit of a chosen path to a window's subgraph, or fall back.
//
// Each reference path is cut into bins of `bin` bp. The nodes of a bin are
// those the reader's context expansion reaches from the reference nodes that
// overlap the bin, with a context of at most `context` bp: the same Dijkstra
// over node sides, from every such node at distance 0, where the reader starts
// from the nodes of its window at distances of 0 or more. So a window's
// subgraph, before its snarls, holds only nodes of the bins the window
// touches. HaplotypeBinNodes lists each bin's nodes, and the reader checks its
// subgraph against the bins it touches.
//
// A path is cut into sections at its visits to the sample's anchor nodes. The
// reader walks a chosen path's section between visits to adjacent anchors of
// its reference path when the span of the two anchors meets the bins it
// touches, from `bound` bp before the first visit to `bound` bp past the
// second, and walks `bound` bp to each side of a visit to an anchor that lies
// up to `bound` bp outside those bins. A visit to a node of a bin is a stray
// for that bin when those walks can miss it. HaplotypeStrays lists the strays
// of each path and bin as stretches of the path, each with the GBWT position
// of its first visit, and the reader walks the rows of chosen paths in the
// bins it touches.
//
// A snarl the reader fills between two boundary nodes of its subgraph adds
// nodes of no bin. A path that passes one of them enters through a boundary
// node, whose visit a walk passes, and the piece through that visit runs on
// into the snarl. A path that lies inside the snarl from end to end has a row
// naming the snarl, listed in the bins of the snarl's lower boundary node.

use super::{fit_u32, Anchors, PathSource};
use crate::encode::encode_ids;
use gbz::{support, Orientation, Pos, ENDMARKER};

use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, BinaryHeap, HashMap, HashSet};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

// A snarl with more nodes between its boundaries is left out, and the reader
// falls back when it fills one.
pub const SNARL_NODES: usize = 1 << 20;
const EMPTY: u64 = 0;

pub const RULE: &str = "HaplotypeBins lists, per reference path of each reference sample with anchors and per bin of haplotype_index_stray_bin bp along it, the nodes within haplotype_index_stray_context bp of a reference node overlapping the bin, by the distance of the context expansion, as runs of node ids (LEB128 gap from the previous run's last id, then run length less one) in parts that each restart from 0, and beside them its stray rows: per bin and path, the stretches of the path holding visits to the bin's nodes that the walks of the keep route can miss: a visit is reached when it lies in a section between visits to adjacent anchors of the bin's reference path whose anchor span meets the bin, or within haplotype_index_stray_bound bp along the path of a visit to such a section's anchor that lies at most that far outside the bin, and a node that is the anchor of two multiples ends no section; consecutive strays more than haplotype_index_stray_gap bp apart start a new row. A row with snarl_high > 0 names a path whose nodes all lie between those two boundary nodes of a top-level snarl, listed in the bins of the lower one";

#[derive(Clone, Copy)]
pub struct Options {
    pub context: usize,
    pub bin: usize,
    pub bound: usize,
    pub gap: usize,
    // The most bytes of one node-list part, so a bin's row fits a b-tree cell.
    pub part_bytes: usize,
}

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Stray {
    pub reference_handle: u32,
    pub bin: u32,
    pub path_handle: u32,
    pub path_start: u32,
    pub snarl_low: u32,
    pub snarl_high: u32,
    pub path_end: u32,
    pub node_handle: u32,
    pub node_offset: u32,
}

#[derive(Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct BinNodes {
    pub reference_handle: u32,
    pub bin: u32,
    pub part: u32,
    pub nodes: Vec<u8>,
}

pub struct Output {
    pub rows: Vec<Stray>,
    pub bins: Vec<BinNodes>,
    pub samples: Vec<String>,
    pub snarls_modeled: bool,
    pub chain_links: usize,
}

impl Output {
    pub fn empty() -> Self {
        Output {
            rows: Vec::new(),
            bins: Vec::new(),
            samples: Vec::new(),
            snarls_modeled: false,
            chain_links: 0,
        }
    }
}

fn pack(reference_handle: usize, bin: usize) -> u64 {
    ((reference_handle as u64 + 1) << 32) | fit_u32("bin", bin) as u64
}

fn unpack(locus: u64) -> (usize, usize) {
    ((locus >> 32) as usize - 1, (locus & 0xffff_ffff) as usize)
}

// Up to two bins per node in place; a node of more spills all of them into a
// shared map.
struct Loci {
    slots: Vec<[AtomicU64; 2]>,
    spilled: Vec<AtomicU64>,
    full: Vec<Mutex<HashMap<usize, Vec<u64>>>>,
}

impl Loci {
    fn new(max_node_id: usize) -> Self {
        Loci {
            slots: (0..=max_node_id)
                .map(|_| [AtomicU64::new(EMPTY), AtomicU64::new(EMPTY)])
                .collect(),
            spilled: (0..=max_node_id / 64).map(|_| AtomicU64::new(0)).collect(),
            full: (0..1024).map(|_| Mutex::new(HashMap::new())).collect(),
        }
    }

    fn is_spilled(&self, id: usize) -> bool {
        self.spilled[id / 64].load(Ordering::Relaxed) & (1 << (id % 64)) != 0
    }

    fn note(&self, id: usize, locus: u64) {
        let mut held = [EMPTY; 2];
        for (i, slot) in self.slots[id].iter().enumerate() {
            loop {
                let current = slot.load(Ordering::Relaxed);
                if current == EMPTY {
                    if slot
                        .compare_exchange(EMPTY, locus, Ordering::Relaxed, Ordering::Relaxed)
                        .is_ok()
                    {
                        return;
                    }
                    continue;
                }
                if current == locus {
                    return;
                }
                held[i] = current;
                break;
            }
        }
        let mut shard = self.full[id % self.full.len()].lock().unwrap();
        let loci = shard.entry(id).or_insert_with(|| held.to_vec());
        if !loci.contains(&locus) {
            loci.push(locus);
        }
        self.spilled[id / 64].fetch_or(1 << (id % 64), Ordering::Relaxed);
    }

    fn of(&self, id: usize) -> Vec<u64> {
        if self.is_spilled(id) {
            return self.full[id % self.full.len()]
                .lock()
                .unwrap()
                .get(&id)
                .cloned()
                .unwrap_or_default();
        }
        self.slots[id]
            .iter()
            .map(|s| s.load(Ordering::Relaxed))
            .filter(|&l| l != EMPTY)
            .collect()
    }
}

struct Visit {
    pos: Pos,
    offset: i64,
    len: i64,
}

fn forward_visits(source: &dyn PathSource, path_handle: usize) -> Vec<Visit> {
    let mut visits = Vec::new();
    let mut pos = source.start(path_handle, Orientation::Forward);
    let mut offset = 0i64;
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (len, next) = source.step(current);
        visits.push(Visit {
            pos: current,
            offset,
            len: len as i64,
        });
        offset += len as i64;
        pos = next;
    }
    visits
}

// The regions the reader's fill of a top-level snarl adds: the nodes reached
// from the first boundary's successors and the second boundary's predecessors
// without crossing either boundary node.
pub struct Snarls {
    // node id -> region + 1, or 0
    region_of: Vec<u32>,
    // region -> its boundary node ids, lower first
    bounds: Vec<(u32, u32)>,
    pub modeled: bool,
    pub links: usize,
}

impl Snarls {
    pub fn none() -> Self {
        Snarls {
            region_of: Vec::new(),
            bounds: Vec::new(),
            modeled: false,
            links: 0,
        }
    }

    // `links` holds each link from both of its ends; `count` is the graph's
    // chain_links tag, which the reader compares with its own.
    pub fn find(source: &dyn PathSource, links: &[(usize, usize)], count: usize) -> Self {
        let started = Instant::now();
        let mut region_of = vec![0u32; source.max_node_id() + 1];
        let mut bounds = Vec::new();
        let mut seen = HashSet::new();
        let mut modeled = true;
        let mut oversized = 0;
        for &(handle, next) in links {
            let key = std::cmp::min(
                (handle, next),
                (support::flip_node(next), support::flip_node(handle)),
            );
            let (a, b) = (support::node_id(handle), support::node_id(next));
            if a == b || !seen.insert(key) {
                continue;
            }
            bounds.push((
                fit_u32("snarl_low", a.min(b)),
                fit_u32("snarl_high", a.max(b)),
            ));
            let region = bounds.len() as u32;
            let mut members = Vec::new();
            let mut stack = vec![handle, support::flip_node(next)];
            'fill: while let Some(current) = stack.pop() {
                for successor in source.successors(current) {
                    let id = support::node_id(successor);
                    if id == a || id == b || region_of[id] == region {
                        continue;
                    }
                    if region_of[id] != 0 {
                        modeled = false;
                        continue;
                    }
                    region_of[id] = region;
                    members.push(id);
                    if members.len() > SNARL_NODES {
                        break 'fill;
                    }
                    stack.push(successor);
                    stack.push(support::flip_node(successor));
                }
            }
            if members.len() > SNARL_NODES {
                for id in members {
                    region_of[id] = 0;
                }
                oversized += 1;
            }
        }
        if !modeled {
            eprintln!("Warning: two top-level snarls share a node, so the stray rows leave snarls out and a query that fills one identifies every walk");
        }
        eprintln!(
            "Filled {} top-level snarls ({} above {} nodes left out) in {:.0} s",
            bounds.len(),
            oversized,
            SNARL_NODES,
            started.elapsed().as_secs_f64()
        );
        Snarls {
            region_of,
            bounds,
            modeled,
            links: count,
        }
    }

    // Paths whose nodes all lie in one region: (path, region, start, offset of the last visit).
    fn inside_paths(&self, source: &dyn PathSource) -> Vec<(usize, usize, Pos, i64)> {
        let mut found = Vec::new();
        if !self.modeled {
            return found;
        }
        for path_handle in 0..source.path_count() {
            let Some(start) = source.start(path_handle, Orientation::Forward) else {
                continue;
            };
            if start.node == ENDMARKER {
                continue;
            }
            let region = self.region_of[support::node_id(start.node)];
            if region == 0 {
                continue;
            }
            let mut pos = Some(start);
            let mut offset = 0i64;
            let mut last = 0i64;
            let mut inside = true;
            while let Some(current) = pos {
                if current.node == ENDMARKER {
                    break;
                }
                if self.region_of[support::node_id(current.node)] != region {
                    inside = false;
                    break;
                }
                let (len, next) = source.step(current);
                last = offset;
                offset += len as i64;
                pos = next;
            }
            if inside {
                found.push((path_handle, region as usize - 1, start, last));
            }
        }
        found
    }
}

#[derive(Clone, Copy)]
enum Cut {
    // reference path, anchor multiple, anchor offset on the reference path
    One(usize, i64, i64),
    Several,
}

// The stretch between a path's visits to adjacent anchors of one reference
// path: the anchors' reference offsets, low and high, and the path offsets of
// the visits to them.
#[derive(Clone, Copy)]
struct Section {
    reference: usize,
    low: i64,
    low_visit: i64,
    high: i64,
    high_visit: i64,
}

struct Sample {
    reference: BTreeSet<usize>,
    anchor: HashMap<usize, Cut>,
    loci: Loci,
    options: Options,
}

impl Sample {
    // The reference nodes overlapping each bin of each reference path.
    fn bins(&self, source: &dyn PathSource) -> Vec<(usize, usize, Vec<usize>)> {
        let mut bins: BTreeMap<(usize, usize), Vec<usize>> = BTreeMap::new();
        for &handle in &self.reference {
            for visit in forward_visits(source, handle) {
                let first = visit.offset as usize / self.options.bin;
                let last = (visit.offset + visit.len.max(1) - 1) as usize / self.options.bin;
                for bin in first..=last {
                    bins.entry((handle, bin))
                        .or_default()
                        .push(support::node_id(visit.pos.node));
                }
            }
        }
        bins.into_iter()
            .map(|((handle, bin), ids)| (handle, bin, ids))
            .collect()
    }

    // Dijkstra over node sides from every source node, with the context
    // expansion's costs: 1 per edge, the node length less one to cross a node.
    fn ball(&self, source: &dyn PathSource, sources: &[usize]) -> Vec<usize> {
        let context = self.options.context as i64;
        let mut heap: BinaryHeap<Reverse<(i64, usize, bool)>> = BinaryHeap::new();
        for &id in sources {
            heap.push(Reverse((0, id, false)));
            heap.push(Reverse((0, id, true)));
        }
        let mut seen: HashSet<(usize, bool)> = HashSet::new();
        let mut ids = BTreeSet::new();
        while let Some(Reverse((distance, id, right))) = heap.pop() {
            if !seen.insert((id, right)) {
                continue;
            }
            ids.insert(id);
            let len = source.node_len(support::encode_node(id, Orientation::Forward)) as i64;
            if !seen.contains(&(id, !right)) && distance + len - 1 <= context {
                heap.push(Reverse((distance + len - 1, id, !right)));
            }
            if distance + 1 <= context {
                let exit = support::encode_node(
                    id,
                    if right {
                        Orientation::Forward
                    } else {
                        Orientation::Reverse
                    },
                );
                for successor in source.successors(exit) {
                    let entry_right = support::node_orientation(successor) == Orientation::Reverse;
                    let successor_id = support::node_id(successor);
                    if !seen.contains(&(successor_id, entry_right)) {
                        heap.push(Reverse((distance + 1, successor_id, entry_right)));
                    }
                }
            }
        }
        ids.into_iter().collect()
    }

    fn section(&self, visits: &[Visit], from: (usize, Cut), to: (usize, Cut)) -> Option<Section> {
        match (from.1, to.1) {
            (Cut::One(path_a, k_a, at_a), Cut::One(path_b, k_b, at_b))
                if path_a == path_b && (k_a - k_b).abs() == 1 =>
            {
                let (visit_a, visit_b) = (visits[from.0].offset, visits[to.0].offset);
                Some(if at_a <= at_b {
                    Section {
                        reference: path_a,
                        low: at_a,
                        low_visit: visit_a,
                        high: at_b,
                        high_visit: visit_b,
                    }
                } else {
                    Section {
                        reference: path_a,
                        low: at_b,
                        low_visit: visit_b,
                        high: at_a,
                        high_visit: visit_a,
                    }
                })
            }
            _ => None,
        }
    }

    fn reached(&self, section: &Section, offset: i64, reference: usize, bin: usize) -> bool {
        let bound = self.options.bound as i64;
        let start = (bin * self.options.bin) as i64;
        let end = start + self.options.bin as i64;
        section.reference == reference
            && ((section.low < end && section.high >= start)
                || (section.low >= end
                    && section.low - end <= bound
                    && (offset - section.low_visit).abs() <= bound)
                || (section.high < start
                    && start - section.high <= bound
                    && (offset - section.high_visit).abs() <= bound))
    }

    fn strays_of(&self, source: &dyn PathSource, path_handle: usize) -> Vec<Stray> {
        let visits = forward_visits(source, path_handle);
        let cuts: Vec<(usize, Cut)> = visits
            .iter()
            .enumerate()
            .filter_map(|(i, visit)| {
                self.anchor
                    .get(&support::node_id(visit.pos.node))
                    .map(|&cut| (i, cut))
            })
            .collect();
        let sections: Vec<Option<Section>> = cuts
            .windows(2)
            .map(|pair| self.section(&visits, pair[0], pair[1]))
            .collect();
        // (bin, visit index)
        let mut strays: Vec<(u64, usize)> = Vec::new();
        let mut before = 0;
        for (i, visit) in visits.iter().enumerate() {
            while before < cuts.len() && cuts[before].0 < i {
                before += 1;
            }
            let at_cut = before < cuts.len() && cuts[before].0 == i;
            let around = [
                if before >= 1 {
                    sections.get(before - 1).copied().flatten()
                } else {
                    None
                },
                if at_cut {
                    sections.get(before).copied().flatten()
                } else {
                    None
                },
            ];
            for locus in self.loci.of(support::node_id(visit.pos.node)) {
                let (reference, bin) = unpack(locus);
                if !around
                    .iter()
                    .flatten()
                    .any(|section| self.reached(section, visit.offset, reference, bin))
                {
                    strays.push((locus, i));
                }
            }
        }
        strays.sort_unstable();
        let mut rows = Vec::new();
        let mut first = 0;
        for next in 1..=strays.len() {
            let split = next == strays.len() || strays[next].0 != strays[first].0 || {
                let previous = &visits[strays[next - 1].1];
                visits[strays[next].1].offset - (previous.offset + previous.len)
                    > self.options.gap as i64
            };
            if split {
                let (reference, bin) = unpack(strays[first].0);
                let start = &visits[strays[first].1];
                rows.push(Stray {
                    reference_handle: reference as u32,
                    bin: bin as u32,
                    path_handle: path_handle as u32,
                    path_start: fit_u32("path_start", start.offset),
                    snarl_low: 0,
                    snarl_high: 0,
                    path_end: fit_u32("path_end", visits[strays[next - 1].1].offset),
                    node_handle: fit_u32("node_handle", start.pos.node),
                    node_offset: fit_u32("node_offset", start.pos.offset),
                });
                first = next;
            }
        }
        rows
    }
}

pub fn reference_samples(
    source: &dyn PathSource,
    anchors: &Anchors,
) -> Vec<(String, BTreeSet<usize>)> {
    let mut by_sample: HashMap<String, BTreeSet<usize>> = HashMap::new();
    for &handle in &anchors.reference_paths {
        by_sample
            .entry(source.sample_of(handle))
            .or_default()
            .insert(handle);
    }
    let mut samples: Vec<(String, BTreeSet<usize>)> = by_sample.into_iter().collect();
    samples.sort();
    samples
}

// Runs `work` over 0..count, on one thread or several sharing one source.
pub trait Runner {
    fn source(&self) -> &dyn PathSource;
    fn run(&self, count: usize, work: &(dyn Fn(&dyn PathSource, usize) + Sync));
}

pub struct Serial<'a>(pub &'a dyn PathSource);

impl Runner for Serial<'_> {
    fn source(&self) -> &dyn PathSource {
        self.0
    }

    fn run(&self, count: usize, work: &(dyn Fn(&dyn PathSource, usize) + Sync)) {
        for i in 0..count {
            work(self.0, i);
        }
    }
}

pub struct Threads<'a> {
    pub source: &'a (dyn PathSource + Sync),
    pub threads: usize,
}

impl Runner for Threads<'_> {
    fn source(&self) -> &dyn PathSource {
        self.source
    }

    fn run(&self, count: usize, work: &(dyn Fn(&dyn PathSource, usize) + Sync)) {
        let next = AtomicUsize::new(0);
        thread::scope(|scope| {
            for _ in 0..self.threads.max(1) {
                scope.spawn(|| loop {
                    let i = next.fetch_add(1, Ordering::Relaxed);
                    if i >= count {
                        break;
                    }
                    work(self.source, i);
                });
            }
        });
    }
}

pub fn strays(
    runner: &dyn Runner,
    anchors: &Anchors,
    spacing: usize,
    options: Options,
    snarls: &Snarls,
) -> Output {
    let mut out = Output::empty();
    out.snarls_modeled = snarls.modeled;
    out.chain_links = snarls.links;
    if options.context == 0 || spacing == 0 {
        return out;
    }
    let source = runner.source();
    let paths = source.path_count();
    let inside = snarls.inside_paths(source);
    for (name, reference) in reference_samples(source, anchors) {
        let started = Instant::now();
        let mut sample = Sample {
            anchor: HashMap::new(),
            loci: Loci::new(source.max_node_id()),
            options,
            reference,
        };
        for a in anchors
            .rows
            .iter()
            .filter(|a| sample.reference.contains(&(a.path_handle as usize)))
        {
            let cut = Cut::One(
                a.path_handle as usize,
                (a.anchor_offset as usize / spacing) as i64,
                a.path_offset as i64,
            );
            sample
                .anchor
                .entry(support::node_id(a.node_handle as usize))
                .and_modify(|held| {
                    if !matches!(*held, Cut::One(path, k, _) if path == a.path_handle as usize && k == (a.anchor_offset as usize / spacing) as i64) {
                        *held = Cut::Several;
                    }
                })
                .or_insert(cut);
        }
        let sample = sample;
        let bins = sample.bins(source);
        let listed = Mutex::new(Vec::new());
        runner.run(bins.len(), &|source, i| {
            let (handle, bin, ids) = &bins[i];
            let nodes = sample.ball(source, ids);
            let locus = pack(*handle, *bin);
            for &id in &nodes {
                sample.loci.note(id, locus);
            }
            let parts: Vec<BinNodes> = encode_ids(&nodes, sample.options.part_bytes)
                .into_iter()
                .enumerate()
                .map(|(part, nodes)| BinNodes {
                    reference_handle: *handle as u32,
                    bin: *bin as u32,
                    part: part as u32,
                    nodes,
                })
                .collect();
            listed.lock().unwrap().extend(parts);
        });
        let listed = listed.into_inner().unwrap();
        eprintln!(
            "Listed the nodes of {} bins for reference sample {} in {:.0} s",
            bins.len(),
            name,
            started.elapsed().as_secs_f64()
        );
        let found = Mutex::new(Vec::new());
        runner.run(paths, &|source, i| {
            let rows = sample.strays_of(source, i);
            if !rows.is_empty() {
                found.lock().unwrap().extend(rows);
            }
        });
        let mut rows = found.into_inner().unwrap();
        let mut snarl_rows = 0;
        for &(path_handle, region, start, last) in &inside {
            let (low, high) = snarls.bounds[region];
            for locus in sample.loci.of(low as usize) {
                let (reference, bin) = unpack(locus);
                rows.push(Stray {
                    reference_handle: reference as u32,
                    bin: bin as u32,
                    path_handle: path_handle as u32,
                    path_start: 0,
                    snarl_low: low,
                    snarl_high: high,
                    path_end: fit_u32("path_end", last),
                    node_handle: fit_u32("node_handle", start.node),
                    node_offset: fit_u32("node_offset", start.offset),
                });
                snarl_rows += 1;
            }
        }
        eprintln!(
            "Found {} stray rows for reference sample {}, {} of them for paths inside a snarl, in {:.0} s",
            rows.len(),
            name,
            snarl_rows,
            started.elapsed().as_secs_f64()
        );
        out.rows.extend(rows);
        out.bins.extend(listed);
        out.samples.push(name);
    }
    out.rows.sort_unstable();
    out.rows.dedup();
    out.bins.sort_unstable();
    out
}
