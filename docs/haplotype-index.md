# Naming haplotypes

Upstream gbz-base names only the query path and prints every other walk as
`unknown#N`. This package names the walks using a haplotype index, a companion
SQLite file that the Rust program `gbz-haplotype-index` writes beside the graph
database. Build the index once per graph, then pass it to the library or the
command line.

The index lives in its own file so the graph database stays exactly what
`gbz-base construct` wrote, and so an index can name the walks of a database
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

Pass `keep` to a range query, or `--keep SAMPLE[#HAP]` on the command line, to
get only the haplotypes you name. `keepHaplotypes(predicate)` leaves the
reference, the chosen walks and the nodes they visit, and drops everything else.

## How walks get their names

The package has two ways to find a walk's name, and the query picks one itself.
Both return the same walks, which `test/anchors.test.ts` checks, so the choice
changes only how long a query takes.

**Sampled**, the default. The query extracts every haplotype through the window,
and `identifyPaths()` names each walk from the samples on the window's nodes. A
walk with no sample in the window is followed past the window for up to four
sampling intervals. If it meets no sample there, it stays unnamed. With `keep`,
the query drops the other walks only after reading and naming all of them, so
the cost follows the number of haplotypes in the graph, 464 in HPRC.

**Anchored**, when a query has `keep` and the index was built with
`--anchor-spacing`. The query reads which haplotypes pass the anchor just before
the window, then walks only the chosen ones forward through it. The cost follows
the number you keep. Keeping 8 HPRC haplotypes takes 0.4-1 s with the pages
cached, against 2-9 s by the sampled route
([performance.md](performance.md#keeping-a-set-of-haplotypes)).

An anchored walk can fail to reach the window, for example when its contig ends
between the anchor and the window. The query then redoes the window by the
sampled route, so the answer stays complete but takes longer.
`gbz-base-query --stats` prints the reason.

`getSubgraphForRange()` and `getAlignmentsForRange()` do all of this for you.
The lower-level queries return unnamed walks until you call `identifyPaths()`.
On the command line, `--resolve` names the walks, and `--alignments`, `--keep`,
`--against` and `--stack` turn it on.

The tests check every name by walking its fragment backward through the GBWT to
the path's recorded start.
