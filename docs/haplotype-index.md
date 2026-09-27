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

For the 10 GB HPRC v2.1 GRCh38 database, the haplotype index is 7.9 GB. Building
it from the GBZ takes 13 minutes on 24 cores and peaks at 12 GB of memory. The
index records the graph's path and node counts, and `open` checks them against
the graph.

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

A pangenome graph like HPRC's holds hundreds of haplotypes, and a query returns
every one that passes through the window, 464 in each HPRC window we measured.
Often you only care about a few, such as the two haplotypes of one sample, or a
handful you want to compare against the reference.

`keep` handles that case. In the library it takes a predicate over each walk's
name; on the command line, `--keep SAMPLE` or `--keep SAMPLE#HAP` picks a sample
or one of its haplotypes, and can repeat. The query then returns the reference,
the chosen haplotypes and only the nodes they visit, so the result is smaller
and, on the anchored route below, much faster to produce.

```bash
gbz-base-query graph.gbz.db --haplotype-index graph.haplotype-index.db \
  --keep HG002 --keep HG00733#1 ...
```

`keep` needs the haplotype index, because the graph database alone cannot tell
which walk belongs to which sample.

## How walks get their names

The GBWT stores each haplotype as an unlabeled walk, and the haplotype index
supplies the names. A query uses it by one of two routes, picking the route
itself. Both return the same walks with the same names and differ only in speed.
The tests check this, and check each name by walking back through the GBWT to
the path's recorded start.

![How a query picks between the sampled and anchored routes](img/naming-routes.svg)

### What a sample is

In the GBWT, each step of a walk is a position: a node, plus an offset saying
which of the walks through that node this one is. From any position the GBWT
gives the next one, so a walk can be followed step by step, but no position says
which path it belongs to.

A sample is one position with that answer written down. Each row of
`HaplotypeSamples` maps a position to its path, the path's orientation and the
path's coordinate there. A row for every step of every haplotype would make the
index far larger than the graph, so `gbz-haplotype-index` writes one every
`--interval` bp along each path (16 kb in the build above), in both
orientations.

To name a walk, a query follows it until it lands on a sampled position. The
walk belongs to that sample's path, and its coordinate is the sample's, shifted
by the bp walked to reach it. A larger interval makes the index smaller and the
walks to the nearest sample longer. The idea is the same as the sampled suffix
array in an FM-index.

An anchor is a sample chosen for coverage instead of spacing. Every
`--anchor-spacing` bp along the reference, `gbz-haplotype-index` picks a node
that most haplotypes pass, and writes a sample for every haplotype's visit
there. Reading the samples at an anchor node lists every haplotype passing that
point, already named.

![The sampled route names every walk in the window from samples; the anchored route walks only the kept haplotypes from the anchor before the window](img/naming-routes-layout.svg)

[naming-routes.dot](img/naming-routes.dot) is the source of the diagram above
the routes, and this one is hand-written SVG.

### Sampled: extract every walk, then name each one

A query without `keep` always takes this route, and so does a `keep` query on an
index without anchors.

1. The query extracts every walk crossing the window's nodes, all of them
   unlabeled.
2. `identifyPaths()` names each walk from a sample on its steps inside the
   window.
3. A walk with no sample in the window is followed past it, up to four sampling
   intervals, and stays unnamed if it meets none.
4. With `keep`, the query drops the unwanted haplotypes only after naming them
   all.

The cost grows with every haplotype in the window, however few `keep` asks for.

### Anchored: start from the names, walk only the kept ones

A `keep` query takes this route when the index was built with
`--anchor-spacing`.

1. The query finds the anchor at or before the window's start, up to one spacing
   (131 kb) upstream, and walks the reference from there to the window.
2. It reads the anchor's rows, which list every haplotype there with its name
   and GBWT step.
3. It keeps the rows for the chosen haplotypes and follows each forward through
   the window. The name is known from the start, so nothing needs identifying.

The other haplotypes are never touched, so the cost grows with the number kept.
With the pages cached, 8 HPRC haplotypes take 0.4-1 s this way, against 2-9 s
sampled ([performance.md](performance.md#keeping-a-set-of-haplotypes)).

### When the anchored route needs help

A kept haplotype can be missing from the anchor's rows, because its contig
starts inside the window or it bypasses the anchor node. The query finds these
from the samples on the window's reference nodes, walks back from each to where
the haplotype enters, and walks forward from there.

A walk from the anchor can also fail to reach the window, for example when its
contig ends first. The anchored result would be incomplete, so the query redoes
the window by the sampled route, and `gbz-base-query --stats` prints which walk
failed and why.
