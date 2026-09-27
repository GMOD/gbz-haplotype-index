# Naming haplotypes

Upstream gbz-base names only the query path and prints every other walk as
`unknown#N`. This package names the walks using a haplotype index, a companion
SQLite file that the Rust program `gbz-haplotype-index` writes beside the graph
database. Build the index once per graph, then pass it to the library or the
command line.

The index is a separate file, so the graph database stays exactly what
`gbz-base construct` wrote, and an index can name the walks of a database
someone else hosts. Its tables could later move into the gbz-base format itself.

## Building the index

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 --anchor-spacing 131072 graph.gbz graph.haplotype-index.db

# without the GBZ, reading paths from the database
gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --from-db graph.gbz.db graph.haplotype-index.db
```

Give `graph.gbz.db` before the index path to check that it matches the GBZ.
`--overwrite` replaces an existing index. The source is in
`tools/haplotype-index/`.

In the library, open the graph and the index together:

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
```

On the command line, pass `--haplotype-index`:

```bash
gbz-base-query https://host/graph.gbz.db \
  --haplotype-index https://host/graph.haplotype-index.db ...
```

For the 10 GB HPRC v2.1 GRCh38 database, the companion is 7.9 GB. Building it
from the GBZ takes 13 minutes on 24 cores and peaks at 12 GB of memory. The
companion records the graph's path and node counts, and `open` checks them
against the graph.

## Tables

- `HaplotypeSamples`: a GBWT position every `--interval` bp along every path, in
  both orientations, with the path and its coordinate there.
- `HaplotypeLengths`: the length of each path.
- `HaplotypeAnchors`: one node every `--anchor-spacing` bp along each reference
  path. `gbz-haplotype-index` picks the node most haplotypes pass through near
  that offset, and stores every haplotype's visit to it as a sample, so the
  samples at an anchor list every haplotype passing that point.

`open` rejects an index built with `--forward-only` with
`ForwardOnlyIndexError`. About half the contigs in a graph like HPRC's are
stored reversed relative to the reference, and walks on those need
reverse-orientation samples.

## Keeping a set of haplotypes

`keep` in the library, or `--keep SAMPLE[#HAP]` on the command line, reduces a
query to the reference, the chosen haplotypes and the nodes they visit.

## How walks get their names

A query names its walks by one of two routes and picks the route itself. The
routes return the same walks and differ only in speed. The tests check this, and
check each name by walking back through the GBWT to the path's recorded start.

**Sampled**, the default. The query extracts every haplotype in the window, and
`identifyPaths()` names each walk from the samples on the window's nodes. A walk
with no sample there is followed up to four sampling intervals past the window,
and stays unnamed if it meets none. With `keep`, the query discards the other
haplotypes afterwards, so the cost grows with every haplotype in the graph.

**Anchored**, for a `keep` query on an index built with `--anchor-spacing`. The
query reads which haplotypes pass the anchor before the window and walks only
the chosen ones, so the cost grows with the number kept. With the pages cached,
8 HPRC haplotypes take 0.4-1 s, against 2-9 s sampled
([performance.md](performance.md#keeping-a-set-of-haplotypes)). If a walk fails
to reach the window, for example because its contig ends first, the query redoes
the window by the sampled route, and `gbz-base-query --stats` prints the reason.
