use gbz::bwt::BWT;
use gbz::support;
use gbz::{GBWT, GBZ, Orientation, Pos, ENDMARKER};
use gbz_base::{GBZBase, GraphInterface};
use rusqlite::{params, Connection};
use simple_sds::serialize;

use std::collections::BTreeSet;
use std::env;
use std::process;
use std::thread;
use std::time::Instant;

const USAGE: &str = "Usage: gbz-haplotype-index [options] graph.gbz graph.gbz.db
       gbz-haplotype-index [options] --output index.db graph.gbz [graph.gbz.db]
       gbz-haplotype-index [options] --from-db graph.gbz.db

Walks every path in both orientations and writes a sample every --interval bp
(default 4096) into table HaplotypeSamples, plus the path lengths into
HaplotypeLengths. The path start and end are always sampled.

Anchors: along every reference path (the paths gbz-base indexes for random
access), one anchor node is chosen per multiple k of --anchor-spacing: the
path's first node for k = 0, and for k >= 1 the node with the most GBWT
positions among those overlapping the half spacing before k * spacing, so it is
one most haplotypes of the region visit. Table HaplotypeAnchors records each
choice, and every path visit through an anchor node is sampled, in both
orientations. A reader can then list every haplotype passing a reference
position, with its own coordinate, from the rows at one node.
--anchor-spacing 0 writes none.

By default the tables are written into the database itself, replacing any
existing ones. With --output FILE they are written into FILE as a standalone
companion database that the reader opens beside the graph database; the
companion records the graph's path and node counts so a mismatch is caught at
open.

With --from-db the walk reads node records from the database itself, so the
GBZ is not needed. Walking a GBZ uses --threads (default: all cores).

Options:
  --interval BP        bp between samples along a path (default 4096)
  --anchor-spacing BP  bp between anchors along a reference path (default 131072)
  --forward-only       sample only the forward orientation of each path
  --output FILE        write a companion database instead of augmenting graph.gbz.db
  --threads N          walker threads for the GBZ route
";

const ANCHOR_RULE: &str = "HaplotypeAnchors names, per indexed reference path and multiple k of the spacing, the path's first node for k = 0 and for k >= 1 the node with the most GBWT positions among those overlapping [k * spacing - spacing / 2, k * spacing), the first on a tie; HaplotypeSamples holds every GBWT position at those nodes, in both orientations";

struct Args {
    gbz: Option<String>,
    db: Option<String>,
    output: Option<String>,
    interval: usize,
    anchor_spacing: usize,
    forward_only: bool,
    threads: usize,
}

fn parse_args() -> Args {
    let mut positional = Vec::new();
    let mut interval = 4096;
    let mut anchor_spacing = 131072;
    let mut forward_only = false;
    let mut from_db = false;
    let mut output = None;
    let mut threads = thread::available_parallelism().map(|n| n.get()).unwrap_or(1);
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
            "--threads" => {
                let value = iter.next().unwrap_or_default();
                threads = value.parse().unwrap_or_else(|_| {
                    eprintln!("Invalid --threads: {}", value);
                    process::exit(1);
                });
            }
            "--output" => output = Some(iter.next().unwrap_or_default()),
            "--forward-only" => forward_only = true,
            "--from-db" => from_db = true,
            "-h" | "--help" => {
                eprint!("{}", USAGE);
                process::exit(0);
            }
            _ => positional.push(arg),
        }
    }
    let valid = if from_db {
        positional.len() == 1
    } else if output.is_some() {
        positional.len() == 1 || positional.len() == 2
    } else {
        positional.len() == 2
    };
    if !valid || interval == 0 || threads == 0 {
        eprint!("{}", USAGE);
        process::exit(1);
    }
    let (gbz, db) = if from_db {
        (None, positional.pop())
    } else {
        let gbz = positional.remove(0);
        (Some(gbz), positional.pop())
    };
    Args { gbz, db, output, interval, anchor_spacing, forward_only, threads }
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

struct NodeSet {
    words: Vec<u64>,
}

impl NodeSet {
    fn with_capacity(max_node_id: usize) -> Self {
        NodeSet { words: vec![0; max_node_id / 64 + 1] }
    }

    fn insert(&mut self, node_id: usize) {
        self.words[node_id / 64] |= 1 << (node_id % 64);
    }

    fn contains(&self, node_id: usize) -> bool {
        self.words.get(node_id / 64).map_or(false, |word| word & (1 << (node_id % 64)) != 0)
    }

    fn len(&self) -> usize {
        self.words.iter().map(|word| word.count_ones() as usize).sum()
    }

    #[cfg(test)]
    fn ids(&self) -> Vec<usize> {
        (0..self.words.len() * 64).filter(|&id| self.contains(id)).collect()
    }
}

trait PathSource {
    fn start(&self, path_handle: usize, orientation: Orientation) -> Option<Pos>;
    fn step(&self, pos: Pos) -> (usize, Option<Pos>);
    fn node_len(&self, handle: usize) -> usize;
    fn visits(&self, handle: usize) -> usize;
    fn indexed_paths(&self) -> Vec<usize>;
    fn max_node_id(&self) -> usize;
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
        bwt.record(index.node_to_record(handle)).map_or(0, |record| record.len())
    }

    fn indexed_paths(&self) -> Vec<usize> {
        let reference_samples: BTreeSet<usize> = self.graph.reference_sample_ids(true).into_iter().collect();
        match self.graph.metadata() {
            Some(metadata) => metadata
                .path_iter()
                .enumerate()
                .filter(|(_, name)| reference_samples.contains(&name.sample()))
                .map(|(handle, _)| handle)
                .collect(),
            None => Vec::new(),
        }
    }

    fn max_node_id(&self) -> usize {
        self.graph.max_node()
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
        Some(if orientation == Orientation::Forward { path.fw_start } else { path.rev_start })
    }

    fn step(&self, pos: Pos) -> (usize, Option<Pos>) {
        let record = self.interface.borrow_mut().get_record(pos.node).unwrap().unwrap();
        (record.sequence_len(), record.to_gbwt_record().lf(pos.offset))
    }

    fn node_len(&self, handle: usize) -> usize {
        self.interface.borrow_mut().get_record(handle).unwrap().unwrap().sequence_len()
    }

    fn visits(&self, handle: usize) -> usize {
        self.interface.borrow_mut().get_record(handle).unwrap().map_or(0, |record| record.to_gbwt_record().len())
    }

    fn indexed_paths(&self) -> Vec<usize> {
        (0..self.paths)
            .filter(|&handle| self.interface.borrow_mut().get_path(handle).unwrap().map_or(false, |path| path.is_indexed))
            .collect()
    }

    fn max_node_id(&self) -> usize {
        self.max_node_id
    }
}

// One anchor per multiple k of the spacing along a reference path: the first
// node for k = 0, and for k >= 1 the node with the most GBWT positions among
// those overlapping [k * spacing - spacing / 2, k * spacing), the first on a
// tie. A window starting at or past k * spacing has that anchor before it, and
// most haplotypes of the region visit it, where the node that happens to
// contain the multiple can be a rare allele.
fn mark_anchors(source: &dyn PathSource, path_handle: usize, spacing: usize, anchors: &mut Vec<Anchor>) {
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
            anchors.push(Anchor { path_handle: path_handle as u32, anchor_offset: 0, node_handle: current.node as u32, path_offset: 0 });
        }
        let end = offset + node_len;
        let mut visits = None;
        loop {
            let target = k * spacing;
            if end <= target - half {
                break;
            }
            let count = *visits.get_or_insert_with(|| source.visits(current.node));
            if best.map_or(true, |(most, _, _)| count > most) {
                best = Some((count, current, offset));
            }
            if end >= target {
                if let Some((_, chosen, chosen_offset)) = best.take() {
                    anchors.push(Anchor { path_handle: path_handle as u32, anchor_offset: target as u32, node_handle: chosen.node as u32, path_offset: chosen_offset as u32 });
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

fn walk(source: &dyn PathSource, path_handle: usize, orientation: Orientation, interval: usize, anchors: &NodeSet, samples: &mut Vec<Sample>) -> usize {
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
                node_handle: current.node as u32,
                node_offset: current.offset as u32,
                path_handle: path_handle as u32,
                orientation: orientation as u8,
                path_offset: offset as u32,
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
            node_handle: end.node as u32,
            node_offset: end.offset as u32,
            path_handle: path_handle as u32,
            orientation: orientation as u8,
            path_offset: end_offset as u32,
        });
    }
    let length = offset;
    if orientation == Orientation::Reverse {
        for sample in samples[first..].iter_mut() {
            let node_len = source.node_len(sample.node_handle as usize);
            sample.path_offset = (length - sample.path_offset as usize - node_len) as u32;
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

fn walk_paths(source: &dyn PathSource, handles: std::ops::Range<usize>, args: &Args, anchors: &NodeSet, label: &str) -> (Vec<Sample>, Vec<(usize, usize)>) {
    let mut samples = Vec::new();
    let mut lengths = Vec::new();
    let started = Instant::now();
    let total = handles.len();
    let mut walked_bp: usize = 0;
    for (done, path_handle) in handles.enumerate() {
        let mut length = 0;
        for &orientation in orientations(args).iter() {
            length = walk(source, path_handle, orientation, args.interval, anchors, &mut samples);
        }
        walked_bp += length;
        lengths.push((path_handle, length));
        if (done + 1) % 500 == 0 || done + 1 == total {
            eprintln!("{}: {} / {} paths, {:.2} Gbp, {} samples, {:.0} s", label, done + 1, total, walked_bp as f64 / 1e9, samples.len(), started.elapsed().as_secs_f64());
        }
    }
    (samples, lengths)
}

struct Anchors {
    nodes: NodeSet,
    rows: Vec<Anchor>,
    reference_paths: usize,
}

fn anchors_gbz(graph: &GBZ, args: &Args) -> Anchors {
    let started = Instant::now();
    let source = GbzSource { graph };
    let reference_paths = source.indexed_paths();
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
    eprintln!("Chose {} anchors on {} distinct nodes over {} reference paths at {} bp spacing in {:.0} s", rows.len(), nodes.len(), reference_paths.len(), args.anchor_spacing, started.elapsed().as_secs_f64());
    Anchors { nodes, rows, reference_paths: reference_paths.len() }
}

fn walk_gbz(graph: &GBZ, paths: usize, args: &Args, anchors: &NodeSet) -> (Vec<Sample>, Vec<(usize, usize)>) {
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
) STRICT;
CREATE TABLE HaplotypeLengths (
    path_handle INTEGER PRIMARY KEY,
    length INTEGER NOT NULL
) STRICT;
CREATE TABLE HaplotypeAnchors (
    path_handle INTEGER NOT NULL,
    anchor_offset INTEGER NOT NULL,
    node_handle INTEGER NOT NULL,
    path_offset INTEGER NOT NULL,
    PRIMARY KEY (path_handle, anchor_offset)
) STRICT;";

fn write(target: &str, standalone: bool, mut samples: Vec<Sample>, lengths: &[(usize, usize)], paths: usize, nodes: usize, anchors: &Anchors, args: &Args) {
    let started = Instant::now();
    samples.sort_unstable_by_key(|s| (s.node_handle, s.node_offset));
    eprintln!("Sorted {} samples in {:.0} s", samples.len(), started.elapsed().as_secs_f64());
    let mut connection = Connection::open(target).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", target, e);
        process::exit(1);
    });
    connection.execute_batch("PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF;").unwrap();
    let mut setup = String::from("DROP TABLE IF EXISTS HaplotypeSamples; DROP TABLE IF EXISTS HaplotypeLengths; DROP TABLE IF EXISTS HaplotypeAnchors; ");
    if standalone {
        setup.push_str("CREATE TABLE IF NOT EXISTS Tags (key TEXT PRIMARY KEY, value TEXT NOT NULL) STRICT; ");
    }
    setup.push_str(SCHEMA);
    connection.execute_batch(&setup).unwrap();
    let transaction = connection.transaction().unwrap();
    {
        let mut write_sample = transaction
            .prepare("INSERT INTO HaplotypeSamples(node_handle, node_offset, path_handle, orientation, path_offset) VALUES (?1, ?2, ?3, ?4, ?5)")
            .unwrap();
        for s in samples.iter() {
            write_sample
                .execute(params![s.node_handle as i64, s.node_offset as i64, s.path_handle as i64, s.orientation as i64, s.path_offset as i64])
                .unwrap();
        }
        let mut write_length = transaction.prepare("INSERT INTO HaplotypeLengths(path_handle, length) VALUES (?1, ?2)").unwrap();
        for &(handle, length) in lengths {
            write_length.execute(params![handle as i64, length as i64]).unwrap();
        }
        let mut write_anchor = transaction
            .prepare("INSERT INTO HaplotypeAnchors(path_handle, anchor_offset, node_handle, path_offset) VALUES (?1, ?2, ?3, ?4)")
            .unwrap();
        let mut rows: Vec<&Anchor> = anchors.rows.iter().collect();
        rows.sort_unstable_by_key(|a| (a.path_handle, a.anchor_offset));
        for a in rows {
            write_anchor.execute(params![a.path_handle as i64, a.anchor_offset as i64, a.node_handle as i64, a.path_offset as i64]).unwrap();
        }
        let mut write_tag = transaction.prepare("INSERT OR REPLACE INTO Tags(key, value) VALUES (?1, ?2)").unwrap();
        let mut drop_tag = transaction.prepare("DELETE FROM Tags WHERE key = ?1").unwrap();
        write_tag.execute(params!["haplotype_index_interval", args.interval.to_string()]).unwrap();
        write_tag.execute(params!["haplotype_index_orientations", if args.forward_only { "forward" } else { "both" }]).unwrap();
        write_tag.execute(params!["haplotype_index_paths", paths.to_string()]).unwrap();
        write_tag.execute(params!["haplotype_index_nodes", nodes.to_string()]).unwrap();
        let anchor_tags = ["haplotype_index_anchor_spacing", "haplotype_index_anchor_rule", "haplotype_index_anchor_paths", "haplotype_index_anchor_nodes"];
        if args.anchor_spacing > 0 {
            let values = [args.anchor_spacing.to_string(), ANCHOR_RULE.to_string(), anchors.reference_paths.to_string(), anchors.nodes.len().to_string()];
            for (key, value) in anchor_tags.iter().zip(values.iter()) {
                write_tag.execute(params![key, value]).unwrap();
            }
        } else {
            for key in anchor_tags.iter() {
                drop_tag.execute(params![key]).unwrap();
            }
        }
    }
    transaction.commit().unwrap();
    eprintln!("Wrote {} samples and {} anchors for {} paths to {} in {:.0} s", samples.len(), anchors.rows.len(), paths, target, started.elapsed().as_secs_f64());
}

fn count_from_db(db: &str, key: &str) -> usize {
    let connection = Connection::open(db).unwrap();
    let value: String = connection
        .query_row("SELECT value FROM Tags WHERE key = ?1", params![key], |row| row.get(0))
        .unwrap_or_else(|_| "0".to_string());
    value.parse().unwrap_or(0)
}

fn max_node_id_from_db(db: &str) -> usize {
    let connection = Connection::open(db).unwrap();
    let handle: i64 = connection.query_row("SELECT max(handle) FROM Nodes", [], |row| row.get(0)).unwrap_or(0);
    support::node_id(handle as usize)
}

fn walk_db(db: &str, args: &Args) -> (Vec<Sample>, Vec<(usize, usize)>, Anchors) {
    let paths = count_from_db(db, "paths");
    let max_node_id = max_node_id_from_db(db);
    let database = GBZBase::open(db).unwrap_or_else(|e| {
        eprintln!("Cannot open {} as a GBZ-base: {}", db, e);
        process::exit(1);
    });
    let interface = GraphInterface::new(&database).unwrap();
    let source = DbSource { interface: std::cell::RefCell::new(interface), paths, max_node_id };
    let reference_paths = source.indexed_paths();
    let rows = anchor_rows(&source, &reference_paths, args.anchor_spacing);
    let nodes = anchor_set(&source, &rows);
    eprintln!("Chose {} anchors on {} distinct nodes over {} reference paths at {} bp spacing", rows.len(), nodes.len(), reference_paths.len(), args.anchor_spacing);
    let (samples, lengths) = walk_paths(&source, 0..paths, args, &nodes, "database walk");
    (samples, lengths, Anchors { nodes, rows, reference_paths: reference_paths.len() })
}

fn main() {
    let args = parse_args();
    let (samples, lengths, paths, nodes, anchors) = match &args.gbz {
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
            eprintln!("Loaded {} with {} paths in {:.0} s", gbz, paths, started.elapsed().as_secs_f64());
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
            let (samples, lengths) = walk_gbz(&graph, paths, &args, &anchors.nodes);
            (samples, lengths, paths, nodes, anchors)
        }
        None => {
            let db = args.db.as_ref().unwrap();
            let paths = count_from_db(db, "paths");
            let nodes = count_from_db(db, "nodes");
            let (samples, lengths, anchors) = walk_db(db, &args);
            (samples, lengths, paths, nodes, anchors)
        }
    };
    match (&args.output, &args.db) {
        (Some(output), _) => write(output, true, samples, &lengths, paths, nodes, &anchors, &args),
        (None, Some(db)) => write(db, false, samples, &lengths, paths, nodes, &anchors, &args),
        (None, None) => unreachable!(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn split_contig() -> String {
        format!("{}/../../test/data/split-contig.gbz.db", env!("CARGO_MANIFEST_DIR"))
    }

    fn args(interval: usize, anchor_spacing: usize) -> Args {
        Args { gbz: None, db: Some(split_contig()), output: None, interval, anchor_spacing, forward_only: false, threads: 1 }
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
        let source = DbSource { interface: std::cell::RefCell::new(interface), paths, max_node_id };
        f(&source)
    }

    #[test]
    fn every_reference_path_gets_its_first_node_and_the_most_visited_node_before_each_multiple() {
        let spacing = 300;
        let (_, _, anchors) = walk_db(&split_contig(), &args(200, spacing));
        assert_eq!(anchors.reference_paths, 2);
        let checked = with_source(|source| {
            let mut checked = 0;
            for &path_handle in source.indexed_paths().iter() {
                let nodes = forward_nodes(source, path_handle);
                let length: usize = nodes.iter().map(|n| n.2).sum();
                let first = anchors.rows.iter().find(|a| a.path_handle as usize == path_handle && a.anchor_offset == 0).unwrap();
                assert_eq!((first.node_handle as usize, first.path_offset), (nodes[0].0, 0));
                let mut k = 1;
                while k * spacing <= length {
                    let target = k * spacing;
                    let overlapping: Vec<&(usize, usize, usize)> =
                        nodes.iter().filter(|(_, start, len)| *start < target && start + len > target - spacing / 2).collect();
                    let most = overlapping.iter().map(|(handle, _, _)| source.visits(*handle)).max().unwrap();
                    let chosen = anchors.rows.iter().find(|a| a.path_handle as usize == path_handle && a.anchor_offset as usize == target).unwrap();
                    let (handle, start, _) = overlapping.iter().find(|(handle, _, _)| *handle == chosen.node_handle as usize).unwrap();
                    assert_eq!(chosen.path_offset as usize, *start);
                    assert_eq!(source.visits(*handle), most);
                    assert!(overlapping.iter().take_while(|(h, _, _)| h != handle).all(|(h, _, _)| source.visits(*h) < most));
                    checked += 1;
                    k += 1;
                }
            }
            checked
        });
        assert!(checked >= 2);
        let distinct: BTreeSet<usize> = anchors.rows.iter().map(|a| support::node_id(a.node_handle as usize)).collect();
        assert_eq!(anchors.nodes.len(), distinct.len());
    }

    #[test]
    fn every_visit_through_an_anchor_node_is_sampled_in_both_orientations() {
        let (samples, lengths, anchors) = walk_db(&split_contig(), &args(200, 300));
        assert_eq!(lengths.len(), 6);
        let mut rows_by_handle = std::collections::BTreeMap::new();
        for sample in samples.iter() {
            if anchors.nodes.contains(support::node_id(sample.node_handle as usize)) {
                rows_by_handle.entry(sample.node_handle as usize).or_insert_with(Vec::new).push(sample.node_offset as usize);
            }
        }
        let expected_handles: BTreeSet<usize> = anchors.nodes.ids().iter().flat_map(|&id| [2 * id, 2 * id + 1]).collect();
        assert_eq!(rows_by_handle.keys().copied().collect::<BTreeSet<_>>(), expected_handles);
        with_source(|source| {
            for (handle, offsets) in rows_by_handle.iter_mut() {
                offsets.sort_unstable();
                assert_eq!(*offsets, (0..source.visits(*handle)).collect::<Vec<_>>());
            }
        });
        let mut keys: Vec<(u32, u32)> = samples.iter().map(|s| (s.node_handle, s.node_offset)).collect();
        let before = keys.len();
        keys.sort_unstable();
        keys.dedup();
        assert_eq!(keys.len(), before);
    }

    #[test]
    fn a_zero_spacing_marks_nothing_and_leaves_the_per_path_samples_alone() {
        let (with, _, anchors) = walk_db(&split_contig(), &args(200, 0));
        assert_eq!(anchors.nodes.len(), 0);
        assert!(anchors.rows.is_empty());
        assert_eq!(with.len(), 42);
    }
}
