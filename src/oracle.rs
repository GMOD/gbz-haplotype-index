// An oracle for the stray rows: the keep route's walk plan, as the proof in
// docs/haplotype-index.md and src/chosenPaths.ts in gbz-base-js state it,
// replayed over every window of whole bins. Every visit of every path to a
// node the window's bins list must lie on a planned walk or in a stray row of
// those bins.

use super::*;
use std::collections::{BTreeMap, BTreeSet, HashMap};

fn decode_part(part: &[u8], ids: &mut BTreeSet<usize>) {
    let mut values = Vec::new();
    let (mut value, mut shift) = (0usize, 0);
    for &byte in part {
        value |= ((byte & 0x7f) as usize) << shift;
        shift += 7;
        if byte & 0x80 == 0 {
            values.push(value);
            value = 0;
            shift = 0;
        }
    }
    let mut previous = 0;
    for pair in values.chunks(2) {
        let start = previous + pair[0];
        ids.extend(start..=start + pair[1]);
        previous = start + pair[1];
    }
}

struct Visit {
    node: usize,
    offset: i64,
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
            node: support::node_id(current.node),
            offset,
        });
        offset += len as i64;
        pos = next;
    }
    visits
}

fn with_source<T>(db: &str, f: impl FnOnce(&DbSource) -> T) -> T {
    let paths = count_from_db(db, "paths");
    let max_node_id = max_node_id_from_db(db);
    let database = GBZBase::open(db).unwrap();
    let interface = GraphInterface::new(&database).unwrap();
    let source = DbSource {
        interface: std::cell::RefCell::new(interface),
        paths,
        max_node_id,
    };
    f(&source)
}

struct Setting {
    interval: usize,
    spacing: usize,
    bin: usize,
    bound: usize,
    context: usize,
    widest_window: usize,
}

fn args(db: &str, s: &Setting) -> Args {
    Args {
        gbz: None,
        db: Some(db.to_string()),
        output: String::new(),
        overwrite: false,
        interval: s.interval,
        anchor_spacing: s.spacing,
        anchor_sample: None,
        reference_interval: s.interval,
        stray_context: s.context,
        stray_bin: s.bin,
        stray_bound: s.bound,
        stray_gap: 10,
        forward_only: false,
        page_size: 4096,
        overview_bin: 0,
        overview_chunk: 256,
        overview_sv: 50,
        threads: 1,
    }
}

// Returns (visits checked, visits found only in a stray row, snarl rows).
fn replay(db: &str, s: &Setting) -> (usize, usize, usize) {
    let args = args(db, s);
    let (_, _, anchors, strays, _) = walk_db(db, &args);
    let spacing = s.spacing as i64;
    let bin_len = s.bin as i64;
    let bound = s.bound as i64;
    let mut listed: BTreeMap<(usize, usize), BTreeSet<usize>> = BTreeMap::new();
    for b in &strays.bins {
        decode_part(
            &b.nodes,
            listed
                .entry((b.reference_handle as usize, b.bin as usize))
                .or_default(),
        );
    }
    let snarl_rows = strays.rows.iter().filter(|r| r.snarl_high > 0).count();
    let mut checked = 0;
    let mut by_row_only = 0;
    with_source(db, |source| {
        let paths: Vec<Vec<Visit>> = (0..source.path_count())
            .map(|p| forward_visits(source, p))
            .collect();
        for &reference in &anchors.reference_paths {
            let anchor_of: BTreeMap<i64, (usize, i64)> = anchors
                .rows
                .iter()
                .filter(|a| a.path_handle as usize == reference)
                .map(|a| {
                    (
                        a.anchor_offset as i64 / spacing,
                        (
                            support::node_id(a.node_handle as usize),
                            a.path_offset as i64,
                        ),
                    )
                })
                .collect();
            let bins: Vec<usize> = listed
                .keys()
                .filter(|(r, _)| *r == reference)
                .map(|(_, b)| *b)
                .collect();
            for &first in &bins {
                for &last in bins
                    .iter()
                    .filter(|&&l| l >= first && l - first < s.widest_window)
                {
                    let lo = first as i64 * bin_len;
                    let hi = (last as i64 + 1) * bin_len;
                    // The reader's anchor reads: chosenPaths.ts, strayRoute.
                    let mut lowest = ((lo - bound).div_euclid(spacing) - 1).max(0);
                    let mut highest = (hi + bound + spacing - 1) / spacing + 1;
                    while lowest > 0 && anchor_of.get(&lowest).is_none_or(|a| a.1 >= lo - bound) {
                        lowest -= 1;
                    }
                    while anchor_of.get(&highest).is_some_and(|a| a.1 <= hi + bound) {
                        highest += 1;
                    }
                    let read: BTreeMap<i64, (usize, i64)> = anchor_of
                        .range(lowest..=highest)
                        .map(|(k, a)| (*k, *a))
                        .collect();
                    if read.is_empty() {
                        continue;
                    }
                    let mut multiples_of: HashMap<usize, Vec<i64>> = HashMap::new();
                    for (k, (node, _)) in &read {
                        multiples_of.entry(*node).or_default().push(*k);
                    }
                    let window_nodes: BTreeSet<usize> = (first..=last)
                        .flat_map(|b| listed.get(&(reference, b)).cloned().unwrap_or_default())
                        .collect();
                    for (path, visits) in paths.iter().enumerate() {
                        let mut walks: Vec<(i64, i64)> = Vec::new();
                        let cuts: Vec<(Option<i64>, i64)> = visits
                            .iter()
                            .filter_map(|v| {
                                multiples_of.get(&v.node).map(|ks| {
                                    (if ks.len() == 1 { Some(ks[0]) } else { None }, v.offset)
                                })
                            })
                            .collect();
                        for pair in cuts.windows(2) {
                            let ((Some(ka), u), (Some(kb), v)) = (pair[0], pair[1]) else {
                                continue;
                            };
                            if (ka - kb).abs() != 1 {
                                continue;
                            }
                            let (a, a_visit, b, b_visit) = if read[&ka].1 <= read[&kb].1 {
                                (read[&ka].1, u, read[&kb].1, v)
                            } else {
                                (read[&kb].1, v, read[&ka].1, u)
                            };
                            if a < hi && b >= lo {
                                walks.push((u - bound, v + bound));
                            } else if a >= hi && a - hi <= bound {
                                walks.push((a_visit - bound, a_visit + bound));
                            } else if b < lo && lo - b <= bound {
                                walks.push((b_visit - bound, b_visit + bound));
                            }
                        }
                        let rows: Vec<(i64, i64)> = strays
                            .rows
                            .iter()
                            .filter(|r| {
                                r.reference_handle as usize == reference
                                    && (first..=last).contains(&(r.bin as usize))
                                    && r.path_handle as usize == path
                                    && r.snarl_high == 0
                            })
                            .map(|r| (r.path_start as i64, r.path_end as i64))
                            .collect();
                        let within = |ranges: &[(i64, i64)], o: i64| {
                            ranges.iter().any(|&(x, y)| x <= o && o <= y)
                        };
                        for v in visits.iter().filter(|v| window_nodes.contains(&v.node)) {
                            let walked = within(&walks, v.offset);
                            assert!(
                                walked || within(&rows, v.offset),
                                "{}: reference {} bins {}..={} path {} visit to node {} at {}: no planned walk {:?} or stray row {:?} reaches it",
                                db, reference, first, last, path, v.node, v.offset, walks, rows
                            );
                            checked += 1;
                            if !walked {
                                by_row_only += 1;
                            }
                        }
                    }
                }
            }
        }
    });
    (checked, by_row_only, snarl_rows)
}

fn fixture(name: &str) -> String {
    format!("{}/test/data/{}.gbz.db", env!("CARGO_MANIFEST_DIR"), name)
}

#[test]
fn every_visit_to_a_window_node_lies_on_a_planned_walk_or_in_a_stray_row() {
    let s = |interval, spacing, bin, bound, context, widest_window| Setting {
        interval,
        spacing,
        bin,
        bound,
        context,
        widest_window,
    };
    let runs = [
        ("split-contig", s(200, 300, 16384, 32768, 1000, 1)),
        ("split-contig", s(200, 300, 100, 150, 50, 4)),
        ("split-contig", s(200, 300, 100, 0, 50, 4)),
        ("split-contig", s(200, 300, 50, 100, 20, 6)),
        ("split-contig", s(200, 300, 1, 100, 1, 3)),
        ("anchor-at-bound", s(1000, 500, 1000, 2000, 1000, 3)),
        ("anchor-at-bound", s(1000, 500, 250, 500, 100, 4)),
        ("keep-bounds", s(4096, 1000, 2000, 1000, 1, 3)),
        ("keep-bounds", s(4096, 1000, 500, 300, 100, 4)),
        ("stray-end", s(256, 2048, 1024, 2048, 100, 3)),
        ("stray-end", s(256, 2048, 256, 512, 100, 4)),
        ("far-pass", s(65536, 16384, 1024, 2048, 100, 3)),
        ("far-pass", s(65536, 4096, 512, 1024, 1000, 3)),
        ("two-copies", s(65536, 16384, 1024, 2048, 100, 3)),
        ("inversion", s(4096, 65536, 1024, 2048, 100, 3)),
        ("inversion", s(4096, 4096, 512, 1024, 1000, 3)),
        ("unplaced", s(65536, 16384, 1024, 2048, 100, 3)),
        ("micb-kir3dl1", s(1000, 2500, 1024, 2048, 100, 3)),
        ("micb-kir3dl1", s(1000, 2500, 300, 500, 1000, 3)),
    ];
    for (name, setting) in &runs {
        let (checked, by_row_only, snarl_rows) = replay(&fixture(name), setting);
        eprintln!(
            "{} spacing {} bin {} bound {} context {}: {} visits checked, {} only in a stray row, {} snarl rows",
            name, setting.spacing, setting.bin, setting.bound, setting.context, checked, by_row_only, snarl_rows
        );
        assert!(checked > 0, "{}", name);
    }
}
