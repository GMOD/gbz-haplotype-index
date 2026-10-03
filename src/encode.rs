// Blob encodings of the index, all sequences of LEB128 varints, cut into parts
// that fit one b-tree cell of the chosen page size without overflow pages. A
// reader concatenates the parts of a key.

use super::Sample;
use crate::strays::Stray;

pub fn push_varint(out: &mut Vec<u8>, mut value: usize) {
    while value >= 0x80 {
        out.push((value & 0x7f) as u8 | 0x80);
        value >>= 7;
    }
    out.push(value as u8);
}

// The most payload a cell of a table b-tree keeps on its page, from the SQLite
// file format, less room for the key columns and the record header. The
// tables with blobs keep their rowid, so their key index holds the keys alone.
pub fn part_limit(page_size: usize) -> usize {
    (page_size - 35).saturating_sub(64)
}

// Sorted distinct ids as runs, in parts of at most `limit` bytes that each
// restart from 0.
pub fn encode_ids(ids: &[usize], limit: usize) -> Vec<Vec<u8>> {
    let mut parts = Vec::new();
    let mut part = Vec::new();
    let mut previous = 0;
    let mut i = 0;
    while i < ids.len() {
        let mut j = i;
        while j + 1 < ids.len() && ids[j + 1] == ids[j] + 1 {
            j += 1;
        }
        let mut run = Vec::new();
        push_varint(&mut run, ids[i] - previous);
        push_varint(&mut run, j - i);
        if !part.is_empty() && part.len() + run.len() > limit {
            parts.push(std::mem::take(&mut part));
            run.clear();
            push_varint(&mut run, ids[i]);
            push_varint(&mut run, j - i);
        }
        part.extend(run);
        previous = ids[j];
        i = j + 1;
    }
    if !part.is_empty() {
        parts.push(part);
    }
    parts
}

fn cut_records(records: Vec<Vec<u8>>, limit: usize) -> Vec<Vec<u8>> {
    let mut parts = Vec::new();
    let mut part = Vec::new();
    for record in records {
        if !part.is_empty() && part.len() + record.len() > limit {
            parts.push(std::mem::take(&mut part));
        }
        part.extend(record);
    }
    if !part.is_empty() {
        parts.push(part);
    }
    parts
}

// The stray rows of one bin: per row, path_handle, path_start, path_end,
// snarl_low, snarl_high, node_handle, node_offset.
pub fn encode_strays(rows: &[Stray], limit: usize) -> Vec<Vec<u8>> {
    let records = rows
        .iter()
        .map(|r| {
            let mut out = Vec::new();
            for value in [
                r.path_handle,
                r.path_start,
                r.path_end,
                r.snarl_low,
                r.snarl_high,
                r.node_handle,
                r.node_offset,
            ] {
                push_varint(&mut out, value as usize);
            }
            out
        })
        .collect();
    cut_records(records, limit)
}

// The forward-orientation samples at both handles of an anchor node: per
// visit, path_handle * 2 + the handle's orientation bit, node_offset and
// path_offset.
pub fn encode_visits(samples: &[Sample], limit: usize) -> Vec<Vec<u8>> {
    let records = samples
        .iter()
        .map(|s| {
            let mut out = Vec::new();
            push_varint(
                &mut out,
                s.path_handle as usize * 2 + (s.node_handle as usize & 1),
            );
            push_varint(&mut out, s.node_offset as usize);
            push_varint(&mut out, s.path_offset as usize);
            out
        })
        .collect();
    cut_records(records, limit)
}

#[cfg(test)]
pub mod tests {
    use super::*;

    pub fn varints(parts: &[Vec<u8>]) -> Vec<usize> {
        let mut values = Vec::new();
        let mut value = 0usize;
        let mut shift = 0;
        for &byte in parts.iter().flatten() {
            value |= ((byte & 0x7f) as usize) << shift;
            shift += 7;
            if byte & 0x80 == 0 {
                values.push(value);
                value = 0;
                shift = 0;
            }
        }
        values
    }

    pub fn decode_ids(parts: &[Vec<u8>]) -> Vec<usize> {
        let mut ids = Vec::new();
        for part in parts {
            let mut previous = 0;
            for pair in varints(std::slice::from_ref(part)).chunks(2) {
                let start = previous + pair[0];
                ids.extend(start..=start + pair[1]);
                previous = start + pair[1];
            }
        }
        ids
    }

    #[test]
    fn ids_round_trip_through_runs_and_parts() {
        let mut ids: Vec<usize> = (5..40).collect();
        ids.extend([100, 101, 5000, 139_000_000, 139_000_001]);
        ids.extend((200_000_000..200_000_000 + 3000).step_by(3));
        ids.sort_unstable();
        let parts = encode_ids(&ids, 800);
        assert!(parts.len() > 1);
        assert!(parts.iter().all(|part| part.len() <= 800));
        assert_eq!(decode_ids(&parts), ids);
        assert!(encode_ids(&[], 800).is_empty());
    }

    #[test]
    fn strays_and_visits_cut_between_records() {
        let rows: Vec<Stray> = (0..50)
            .map(|i| Stray {
                reference_handle: 0,
                bin: 1,
                path_handle: 1000 + i,
                path_start: 3_000_000 * i,
                snarl_low: 0,
                snarl_high: 0,
                path_end: 3_000_000 * i + 500,
                node_handle: 2 * i,
                node_offset: i,
            })
            .collect();
        let parts = encode_strays(&rows, 100);
        assert!(parts.len() > 1);
        assert!(parts.iter().all(|part| part.len() <= 100));
        let values = varints(&parts);
        assert_eq!(values.len(), 7 * rows.len());
        assert_eq!(values[7 * 49], 1049);
        assert_eq!(values[7 * 49 + 2], 3_000_000 * 49 + 500);
        let samples: Vec<Sample> = (0..30)
            .map(|i| Sample {
                node_handle: 1_000_000 + (i % 2),
                node_offset: i,
                path_handle: 40_000 + i,
                orientation: 0,
                path_offset: 100_000_000 + i,
            })
            .collect();
        let parts = encode_visits(&samples, 60);
        assert!(parts.len() > 1);
        let values = varints(&parts);
        assert_eq!(values.len(), 3 * samples.len());
        assert_eq!(values[3], (40_001) * 2 + 1);
        assert_eq!(values[5], 100_000_001);
    }

    #[test]
    fn the_part_limit_follows_the_page_size() {
        assert_eq!(part_limit(4096), 3997);
        assert_eq!(part_limit(65536), 65437);
    }
}
