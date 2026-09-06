use gbz::bwt::BWT;
use gbz::support;
use gbz::{GBWT, GBZ, Orientation, Pos, ENDMARKER};
use gbz_base::{GBZBase, GraphInterface};
use rusqlite::{params, Connection};
use simple_sds::serialize;

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

By default the tables are written into the database itself, replacing any
existing ones. With --output FILE they are written into FILE as a standalone
companion database that the reader opens beside the graph database; the
companion records the graph's path count so a mismatch is caught at open.

With --from-db the walk reads node records from the database itself, so the
GBZ is not needed. Walking a GBZ uses --threads (default: all cores).

Options:
  --interval BP     bp between samples along a path (default 4096)
  --forward-only    sample only the forward orientation of each path
  --output FILE     write a companion database instead of augmenting graph.gbz.db
  --threads N       walker threads for the GBZ route
";

struct Args {
    gbz: Option<String>,
    db: Option<String>,
    output: Option<String>,
    interval: usize,
    forward_only: bool,
    threads: usize,
}

fn parse_args() -> Args {
    let mut positional = Vec::new();
    let mut interval = 4096;
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
    Args { gbz, db, output, interval, forward_only, threads }
}

#[derive(Clone, Copy)]
struct Sample {
    node_handle: u32,
    node_offset: u32,
    path_handle: u32,
    orientation: u8,
    path_offset: u32,
}

trait PathSource {
    fn start(&self, path_handle: usize, orientation: Orientation) -> Option<Pos>;
    fn step(&self, pos: Pos) -> (usize, Option<Pos>);
    fn node_len(&self, handle: usize) -> usize;
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
}

struct DbSource<'a> {
    interface: std::cell::RefCell<GraphInterface<'a>>,
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
}

fn walk(source: &dyn PathSource, path_handle: usize, orientation: Orientation, interval: usize, samples: &mut Vec<Sample>) -> usize {
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
        if offset >= next_sample {
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

fn walk_paths(source: &dyn PathSource, handles: std::ops::Range<usize>, args: &Args, label: &str) -> (Vec<Sample>, Vec<(usize, usize)>) {
    let mut samples = Vec::new();
    let mut lengths = Vec::new();
    let started = Instant::now();
    let total = handles.len();
    let mut walked_bp: usize = 0;
    for (done, path_handle) in handles.enumerate() {
        let mut length = 0;
        for &orientation in orientations(args).iter() {
            length = walk(source, path_handle, orientation, args.interval, &mut samples);
        }
        walked_bp += length;
        lengths.push((path_handle, length));
        if (done + 1) % 500 == 0 || done + 1 == total {
            eprintln!("{}: {} / {} paths, {:.2} Gbp, {} samples, {:.0} s", label, done + 1, total, walked_bp as f64 / 1e9, samples.len(), started.elapsed().as_secs_f64());
        }
    }
    (samples, lengths)
}

fn walk_gbz(graph: &GBZ, paths: usize, args: &Args) -> (Vec<Sample>, Vec<(usize, usize)>) {
    let chunk = (paths + args.threads - 1) / args.threads;
    let started = Instant::now();
    let results: Vec<(Vec<Sample>, Vec<(usize, usize)>)> = thread::scope(|scope| {
        let workers: Vec<_> = (0..args.threads)
            .map(|t| {
                let range = (t * chunk).min(paths)..((t + 1) * chunk).min(paths);
                scope.spawn(move || {
                    let source = GbzSource { graph };
                    let label = format!("thread {} (paths {}..{})", t, range.start, range.end);
                    let result = walk_paths(&source, range.clone(), args, &label);
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
) STRICT;";

fn write(target: &str, standalone: bool, mut samples: Vec<Sample>, lengths: &[(usize, usize)], paths: usize, args: &Args) {
    let started = Instant::now();
    samples.sort_unstable_by_key(|s| (s.node_handle, s.node_offset));
    eprintln!("Sorted {} samples in {:.0} s", samples.len(), started.elapsed().as_secs_f64());
    let mut connection = Connection::open(target).unwrap_or_else(|e| {
        eprintln!("Cannot open {}: {}", target, e);
        process::exit(1);
    });
    connection.execute_batch("PRAGMA journal_mode = OFF; PRAGMA synchronous = OFF;").unwrap();
    let mut setup = String::from("DROP TABLE IF EXISTS HaplotypeSamples; DROP TABLE IF EXISTS HaplotypeLengths; ");
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
        let mut write_tag = transaction.prepare("INSERT OR REPLACE INTO Tags(key, value) VALUES (?1, ?2)").unwrap();
        write_tag.execute(params!["haplotype_index_interval", args.interval.to_string()]).unwrap();
        write_tag.execute(params!["haplotype_index_orientations", if args.forward_only { "forward" } else { "both" }]).unwrap();
        write_tag.execute(params!["haplotype_index_paths", paths.to_string()]).unwrap();
    }
    transaction.commit().unwrap();
    eprintln!("Wrote {} samples for {} paths to {} in {:.0} s", samples.len(), paths, target, started.elapsed().as_secs_f64());
}

fn path_count_from_db(db: &str) -> usize {
    let connection = Connection::open(db).unwrap();
    let value: String = connection
        .query_row("SELECT value FROM Tags WHERE key = 'paths'", [], |row| row.get(0))
        .unwrap_or_else(|_| "0".to_string());
    value.parse().unwrap_or(0)
}

fn main() {
    let args = parse_args();
    let (samples, lengths, paths) = match &args.gbz {
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
            if let Some(db) = &args.db {
                let db_paths = path_count_from_db(db);
                if db_paths != paths {
                    eprintln!("{} has {} paths but {} has {}", gbz, paths, db, db_paths);
                    process::exit(1);
                }
            }
            let (samples, lengths) = walk_gbz(&graph, paths, &args);
            (samples, lengths, paths)
        }
        None => {
            let db = args.db.as_ref().unwrap();
            let paths = path_count_from_db(db);
            let database = GBZBase::open(db).unwrap_or_else(|e| {
                eprintln!("Cannot open {} as a GBZ-base: {}", db, e);
                process::exit(1);
            });
            let interface = GraphInterface::new(&database).unwrap();
            let source = DbSource { interface: std::cell::RefCell::new(interface) };
            let (samples, lengths) = walk_paths(&source, 0..paths, &args, "database walk");
            (samples, lengths, paths)
        }
    };
    match (&args.output, &args.db) {
        (Some(output), _) => write(output, true, samples, &lengths, paths, &args),
        (None, Some(db)) => write(db, false, samples, &lengths, paths, &args),
        (None, None) => unreachable!(),
    }
}
