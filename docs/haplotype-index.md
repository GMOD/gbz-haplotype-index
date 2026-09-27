# The haplotype index

Upstream gbz-base prints each walk other than the query path as `unknown#N`.
This package reports the sample, haplotype and contig of every walk, looked up
in a haplotype index: a sidecar SQLite file that the Rust program
`gbz-haplotype-index` writes beside the graph database. The graph database stays
as `gbz-base construct` wrote it, so a haplotype index also works with a
database someone else hosts.

Two things happen at different times. Building the haplotype index happens once
per graph, and the file covers every haplotype in the graph. Identifying walks
happens on every query that reads the index. A caller who wants a subset of the
haplotypes passes the `keep` option with the query
([below](#querying-a-subset-of-the-haplotypes)); the index needs no rebuild for
a different choice.

## Samples and anchors

A path is one contig of one haplotype, named `sample#haplotype#contig`. The GBWT
stores each path as a walk, a sequence of positions: a node, plus the rank of
this walk among the walks through that node. The GBWT maps each position to the
next one, and identifies a path at the start of its walk, which can be a whole
chromosome away from a query window. The haplotype index records the path at
positions along the way:

- A **sample** records a position, the path through it, the orientation and the
  coordinate along that path. `gbz-haplotype-index` writes one every
  `--interval` bp along each path, in both orientations, plus one at the start
  and end of each path. To identify a walk, the library follows the walk to the
  nearest sample. A larger interval gives a smaller index and longer walks, the
  same trade as the sampled suffix array in an FM-index.
- An **anchor** is a reference node that most haplotypes in the region visit,
  one every `--anchor-spacing` bp along each reference path. The haplotype index
  contains a sample for every visit to an anchor node, so the samples at that
  node list every haplotype passing it.

![graph.gbz.db lists the walks at each node by rank and names each path at the start of its walk. The haplotype index adds a sample every --interval bp that maps a position to a path, so a query names a walk from the next sample along it](img/haplotype-samples.svg)

Table `HaplotypeSamples` contains the samples, `HaplotypeAnchors` lists the
anchor nodes, and `HaplotypeLengths` lists the length of each path.

The word "sample" also means an individual, such as `HG002` in a path name. The
rest of this page uses it for the index entry.

## Building the haplotype index

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 --anchor-spacing 131072 graph.gbz graph.haplotype-index.db

# without the GBZ, reading paths from the database
gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --from-db graph.gbz.db graph.haplotype-index.db
```

Give `graph.gbz.db` before the index path to check that it matches the GBZ, and
`--overwrite` to replace an existing index. `--anchor-spacing 0` writes no
anchors, and every query then takes the sampled route below. The source is in
`tools/haplotype-index/`.

Keep the default of sampling both orientations. `open` rejects an index built
with `--forward-only`, because about half the contigs in a graph like HPRC's run
reversed relative to the reference, and identifying a walk on one of those needs
reverse-orientation samples.

For the 10 GB HPRC v2.1 GRCh38 database, the index is 7.9 GB. Building it from
the GBZ takes 13 minutes on 24 cores and peaks at 12 GB of memory.

## Using the haplotype index

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
```

```bash
gbz-base-query https://host/graph.gbz.db \
  --haplotype-index https://host/graph.haplotype-index.db ...
```

`open` checks the path and node counts recorded in the haplotype index against
the graph.

## Querying a subset of the haplotypes

A query returns every haplotype that passes through the window, 464 in each HPRC
window we measured. To get a few of them, such as the two haplotypes of HG002,
pass `keep` with the query:

```ts
const records = await db.getAlignmentsForRange('GRCh38#0#chr6', start, end, {
  keep: name => name.sample === 'HG002',
})
```

```bash
gbz-base-query ... --keep HG002 --keep HG00733#1
```

The query then returns the reference, the haplotypes you kept and the nodes they
visit. `--keep` takes a sample or `sample#haplotype` and can repeat. A query
that uses the `keep` option needs the haplotype index, because the library looks
up the haplotype of each walk in it.

## How a query identifies walks

The library identifies walks by one of two routes and selects the route from the
options passed with the query. Both routes find the same haplotype for every
walk. The tests compare the two, and check each result by walking back through
the GBWT to the start of the path.

![The sampled route extracts every walk through the window and follows each one to a sample. The anchored route reads the samples at one anchor node, which list every haplotype passing it, then follows the chosen walks through the window](img/naming-routes-layout.svg)

- **Sampled**: the query extracts every walk through the window, then follows
  each walk to a sample to identify it. The cost grows with the number of
  haplotypes in the window.
- **Anchored**: the query reads the samples at the anchor node before the
  window, which list every haplotype passing it, then follows the chosen walks.
  The cost grows with the number of haplotypes chosen.

A query that uses the `keep` option takes the anchored route, and every other
query takes the sampled route. The library sends a query that uses the `keep`
option down the sampled route when its `haplotypes` setting is other than
`'all'`, the default, or when `gbz-haplotype-index` built the haplotype index
with `--anchor-spacing 0`. On that route the library identifies every walk, then
drops the haplotypes the predicate rejects.

![A query that uses the keep option and leaves haplotypes at 'all' takes the anchored route; every other query takes the sampled route](img/naming-routes.svg)

The flowchart source is [naming-routes.dot](img/naming-routes.dot); the
schematic is hand-written SVG.

### Sampled

1. The query extracts every walk through the nodes in the window.
2. `identifyPaths()` identifies each walk from a sample the walk passes inside
   the window. A walk with no sample there is followed past the window, for
   about four sampling intervals, and printed as `unknown#N` if that finds none.
3. For a query that uses the `keep` option, the library then drops the
   haplotypes the predicate rejects.

### Anchored

1. The query looks up the anchor for the last multiple of `--anchor-spacing` at
   or before the window's start. That node lies under 1.5 spacings before the
   window, and the query walks the reference from it to the window.
2. The query reads the samples at the anchor node. Each sample names the
   haplotype passing there, so the query keeps the samples of the chosen
   haplotypes and ignores the rest.
3. Starting from each kept sample, the query follows the walk of that haplotype
   forward through the window. The sample already identifies the walk, so no
   sample search follows.
4. A chosen haplotype has no sample at the anchor when the contig starts after
   the anchor or the walk bypasses that node. The query then looks for a sample
   of that haplotype on the reference nodes in the window, walks back from that
   sample to where the haplotype joins the reference, then forward as in step 3.
5. If a walk never reaches the window, for example because the contig ends
   first, or a haplotype in step 4 has no such sample, the query redoes the
   window by the sampled route. `gbz-base-query --stats` prints the reason.

With the pages cached, 8 HPRC haplotypes take 0.4-1 s by this route, against 2-9
s for the sampled query of all 464
([performance.md](performance.md#a-subset-of-the-haplotypes)).
