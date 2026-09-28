// Table HaplotypeStrays: the visits the keep route's walks from anchor visits
// do not reach, for each reference sample with anchors.
//
// Every node gets loci on the reference sample's paths: the position of each
// reference node within `context` bp of it in the graph, by the distance the
// reader's context expansion uses, and the reference visits on both sides of it
// along every path that passes it. A window's subgraph holds a node only near
// one of its loci. A path is cut into sections at its visits to the sample's
// anchors. A visit is reached when each locus of its node lies inside the
// reference span of its section's two visits, those visits being to adjacent
// anchors of one reference path, or within BOUND of that span on the reference
// and within BOUND of that end's visit along the path. The reader walks every
// such section of a chosen path, BOUND past both ends.
// Every other visit is a stray. Strays cluster into rows of at most CHUNK bp on
// the reference and along the path, each with the GBWT position of its first
// visit, so the reader walks a chosen path's rows forward from there.

use super::{Anchors, PathSource};
use gbz::{support, Orientation, Pos, ENDMARKER};

use std::cmp::Reverse;
use std::collections::{BTreeSet, BinaryHeap, HashMap, HashSet};
use std::sync::atomic::{AtomicU64, AtomicUsize, Ordering};
use std::sync::Mutex;
use std::thread;
use std::time::Instant;

pub const BOUND: i64 = 32768;
pub const CHUNK: i64 = 16384;
const EMPTY: u64 = 0;

pub const RULE: &str = "HaplotypeStrays lists, per reference sample with anchors, every visit of a path to a node with a locus its walks from anchor visits do not reach: a node's loci are the reference positions within the stray context of it in the graph and the reference visits beside it along each path; a visit is reached when each locus lies in the reference span of its section's two visits to adjacent anchors of one reference path, or within 32768 bp of that span and of that end's visit along the path; rows cluster strays within 16384 bp on the reference and along the path, with the GBWT position of the first";

#[derive(Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct Stray {
    pub reference_handle: u32,
    pub reference_start: u32,
    pub reference_end: u32,
    pub path_handle: u32,
    pub path_start: u32,
    pub path_end: u32,
    pub node_handle: u32,
    pub node_offset: u32,
}

fn pack(reference_handle: usize, offset: usize) -> u64 {
    ((reference_handle as u64 + 1) << 32) | offset as u64
}

fn unpack(locus: u64) -> (usize, i64) {
    ((locus >> 32) as usize - 1, (locus & 0xffff_ffff) as i64)
}

fn close(a: u64, b: u64) -> bool {
    let (ha, oa) = unpack(a);
    let (hb, ob) = unpack(b);
    ha == hb && (oa - ob).abs() <= BOUND
}

// Up to two loci per node in place; a node with more spills all of them into
// a shared map. Loci within BOUND of one another on a path count once.
struct Loci {
    slots: Vec<[AtomicU64; 2]>,
    spilled: Vec<AtomicU64>,
    full: Vec<Mutex<HashMap<usize, Vec<u64>>>>,
}

impl Loci {
    fn new(max_node_id: usize) -> Self {
        Loci {
            slots: (0..=max_node_id).map(|_| [AtomicU64::new(EMPTY), AtomicU64::new(EMPTY)]).collect(),
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
                    if slot.compare_exchange(EMPTY, locus, Ordering::Relaxed, Ordering::Relaxed).is_ok() {
                        return;
                    }
                    continue;
                }
                if close(current, locus) {
                    return;
                }
                held[i] = current;
                break;
            }
        }
        let mut shard = self.full[id % self.full.len()].lock().unwrap();
        let loci = shard.entry(id).or_insert_with(|| held.to_vec());
        if !loci.iter().any(|&l| close(l, locus)) {
            loci.push(locus);
        }
        self.spilled[id / 64].fetch_or(1 << (id % 64), Ordering::Relaxed);
    }

    fn of(&self, id: usize) -> Vec<u64> {
        if self.is_spilled(id) {
            return self.full[id % self.full.len()].lock().unwrap().get(&id).cloned().unwrap_or_default();
        }
        self.slots[id].iter().map(|s| s.load(Ordering::Relaxed)).filter(|&l| l != EMPTY).collect()
    }
}

struct Sample {
    reference: BTreeSet<usize>,
    // node id -> locus of the reference visit
    position: HashMap<usize, u64>,
    // node id -> (reference path, anchor multiple, locus)
    anchor: HashMap<usize, (usize, i64, u64)>,
    segments: Vec<Vec<(usize, u64)>>,
    loci: Loci,
    context: i64,
}

fn forward_visits(source: &dyn PathSource, path_handle: usize) -> Vec<(Pos, i64)> {
    let mut visits = Vec::new();
    let mut pos = source.start(path_handle, Orientation::Forward);
    let mut offset = 0i64;
    while let Some(current) = pos {
        if current.node == ENDMARKER {
            break;
        }
        let (len, next) = source.step(current);
        visits.push((current, offset));
        offset += len as i64;
        pos = next;
    }
    visits
}

impl Sample {
    // The reference nodes in each BOUND stretch of each reference path, and
    // the locus of every reference node.
    fn map_reference(&mut self, source: &dyn PathSource) {
        for &handle in &self.reference {
            let mut current: Vec<(usize, u64)> = Vec::new();
            let mut start = 0;
            for (pos, offset) in forward_visits(source, handle) {
                if offset - start >= BOUND && !current.is_empty() {
                    self.segments.push(std::mem::take(&mut current));
                    start = offset;
                }
                let id = support::node_id(pos.node);
                let locus = pack(handle, offset as usize);
                self.position.entry(id).or_insert(locus);
                current.push((id, locus));
            }
            if !current.is_empty() {
                self.segments.push(current);
            }
        }
    }

    // Dijkstra over node sides from every reference node of one stretch, with
    // the context expansion's costs: 1 per edge, the node length less one to
    // cross a node.
    fn graph_loci(&self, source: &dyn PathSource, segment: &[(usize, u64)]) {
        let mut heap: BinaryHeap<Reverse<(i64, usize, bool, u64)>> = BinaryHeap::new();
        for &(id, locus) in segment {
            heap.push(Reverse((0, id, false, locus)));
            heap.push(Reverse((0, id, true, locus)));
        }
        let mut seen: HashSet<(usize, bool)> = HashSet::new();
        while let Some(Reverse((distance, id, right, locus))) = heap.pop() {
            if !seen.insert((id, right)) {
                continue;
            }
            self.loci.note(id, locus);
            let len = source.node_len(support::encode_node(id, Orientation::Forward)) as i64;
            if !seen.contains(&(id, !right)) && distance + len - 1 <= self.context {
                heap.push(Reverse((distance + len - 1, id, !right, locus)));
            }
            if distance + 1 <= self.context {
                let exit = support::encode_node(id, if right { Orientation::Forward } else { Orientation::Reverse });
                for successor in source.successors(exit) {
                    let entry_right = support::node_orientation(successor) == Orientation::Reverse;
                    let successor_id = support::node_id(successor);
                    if !seen.contains(&(successor_id, entry_right)) {
                        heap.push(Reverse((distance + 1, successor_id, entry_right, locus)));
                    }
                }
            }
        }
    }

    // The reference visits on both sides of each visit along one path.
    fn flank_loci(&self, source: &dyn PathSource, path_handle: usize) {
        let visits = forward_visits(source, path_handle);
        let ids: Vec<usize> = visits.iter().map(|(pos, _)| support::node_id(pos.node)).collect();
        let mut last = EMPTY;
        for &id in &ids {
            if let Some(&locus) = self.position.get(&id) {
                last = locus;
            }
            if last != EMPTY {
                self.loci.note(id, last);
            }
        }
        let mut next = EMPTY;
        for &id in ids.iter().rev() {
            if let Some(&locus) = self.position.get(&id) {
                next = locus;
            }
            if next != EMPTY {
                self.loci.note(id, next);
            }
        }
    }

    fn strays_of(&self, source: &dyn PathSource, path_handle: usize) -> Vec<Stray> {
        let visits = forward_visits(source, path_handle);
        let cuts: Vec<(usize, (usize, i64, u64))> = visits
            .iter()
            .enumerate()
            .filter_map(|(i, (pos, _))| self.anchor.get(&support::node_id(pos.node)).map(|&a| (i, a)))
            .collect();
        // (locus, path offset, visit index)
        let mut strays: Vec<(u64, i64, usize)> = Vec::new();
        let mut section = 0;
        for (i, &(pos, offset)) in visits.iter().enumerate() {
            while section < cuts.len() && cuts[section].0 < i {
                section += 1;
            }
            if section < cuts.len() && cuts[section].0 == i {
                continue;
            }
            let before = if section > 0 { Some((visits[cuts[section - 1].0].1, cuts[section - 1].1)) } else { None };
            let after = cuts.get(section).map(|&(j, a)| (visits[j].1, a));
            // (reference path, low end, its visit's offset, high end, its visit's offset)
            let span = match (before, after) {
                (Some((pa, a)), Some((pb, b))) if a.0 == b.0 && (a.1 - b.1).abs() == 1 => {
                    let (oa, ob) = (unpack(a.2).1, unpack(b.2).1);
                    Some(if oa <= ob { (a.0, oa, pa, ob, pb) } else { (a.0, ob, pb, oa, pa) })
                }
                _ => None,
            };
            for locus in self.loci.of(support::node_id(pos.node)) {
                let (handle, at) = unpack(locus);
                let reached = span.map_or(false, |(h, lo, plo, hi, phi)| {
                    h == handle
                        && ((lo <= at && at <= hi)
                            || (at < lo && lo - at <= BOUND && (offset - plo).abs() <= BOUND)
                            || (at > hi && at - hi <= BOUND && (offset - phi).abs() <= BOUND))
                });
                if !reached {
                    strays.push((locus, offset, i));
                }
            }
        }
        rows(path_handle, &visits, strays)
    }
}

// Strays sorted by locus, cut where the reference path changes, the gap exceeds
// BOUND or the extent CHUNK, then each run cut into CHUNK bp of path.
fn rows(path_handle: usize, visits: &[(Pos, i64)], mut strays: Vec<(u64, i64, usize)>) -> Vec<Stray> {
    strays.sort_unstable();
    let mut out = Vec::new();
    let mut start = 0;
    for end in 1..=strays.len() {
        let (h0, o0) = unpack(strays[start].0);
        let split = end == strays.len() || {
            let (h, o) = unpack(strays[end].0);
            let (_, prev) = unpack(strays[end - 1].0);
            h != h0 || o - prev > BOUND || o - o0 > CHUNK
        };
        if !split {
            continue;
        }
        let mut run: Vec<(i64, i64, usize)> = strays[start..end].iter().map(|&(l, p, i)| (p, unpack(l).1, i)).collect();
        run.sort_unstable();
        let mut first = 0;
        for last in 1..=run.len() {
            if last == run.len() || run[last].0 - run[first].0 > CHUNK {
                let part = &run[first..last];
                let (pos, _) = visits[part[0].2];
                out.push(Stray {
                    reference_handle: h0 as u32,
                    reference_start: part.iter().map(|r| r.1).min().unwrap() as u32,
                    reference_end: part.iter().map(|r| r.1).max().unwrap() as u32,
                    path_handle: path_handle as u32,
                    path_start: part[0].0 as u32,
                    path_end: part[part.len() - 1].0 as u32,
                    node_handle: pos.node as u32,
                    node_offset: pos.offset as u32,
                });
                first = last;
            }
        }
        start = end;
    }
    out.sort_unstable();
    out.dedup();
    out
}

fn reference_samples(source: &dyn PathSource, anchors: &Anchors) -> Vec<(String, BTreeSet<usize>)> {
    let mut by_sample: HashMap<String, BTreeSet<usize>> = HashMap::new();
    for &handle in &anchors.reference_paths {
        by_sample.entry(source.sample_of(handle)).or_default().insert(handle);
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

pub fn strays(runner: &dyn Runner, anchors: &Anchors, spacing: usize, context: usize) -> (Vec<Stray>, Vec<String>) {
    let mut out = Vec::new();
    let mut names = Vec::new();
    if context == 0 || spacing == 0 {
        return (out, names);
    }
    let source = runner.source();
    let paths = source.path_count();
    for (name, reference) in reference_samples(source, anchors) {
        let started = Instant::now();
        let mut sample = Sample {
            anchor: HashMap::new(),
            position: HashMap::new(),
            segments: Vec::new(),
            loci: Loci::new(source.max_node_id()),
            context: context as i64,
            reference,
        };
        sample.map_reference(source);
        for a in anchors.rows.iter().filter(|a| sample.reference.contains(&(a.path_handle as usize))) {
            let locus = pack(a.path_handle as usize, a.path_offset as usize);
            sample.anchor.entry(support::node_id(a.node_handle as usize)).or_insert((a.path_handle as usize, a.anchor_offset as i64 / spacing as i64, locus));
        }
        let sample = sample;
        runner.run(sample.segments.len(), &|source, i| sample.graph_loci(source, &sample.segments[i]));
        runner.run(paths, &|source, i| sample.flank_loci(source, i));
        let found = Mutex::new(Vec::new());
        runner.run(paths, &|source, i| {
            let rows = sample.strays_of(source, i);
            if !rows.is_empty() {
                found.lock().unwrap().extend(rows);
            }
        });
        let mut rows = found.into_inner().unwrap();
        rows.sort_unstable();
        eprintln!("Found {} stray rows for reference sample {} in {:.0} s", rows.len(), name, started.elapsed().as_secs_f64());
        out.extend(rows);
        names.push(name);
    }
    (out, names)
}
