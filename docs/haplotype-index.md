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
stores each path as a walk, a sequence of positions. A position is a node plus
the rank of one walk among the walks through that node, so each position belongs
to exactly one walk. The GBWT maps each position to the next one, and identifies
a path at the start of its walk, which can be a whole chromosome away from a
query window. The haplotype index records the path at positions along the way:

- A **sample** records a position, the path through it, the orientation and the
  coordinate along that path. `gbz-haplotype-index` writes one every
  `--interval` bp along each path, in both orientations, plus one at the start
  and end of each path. To identify a walk, the library follows the walk to the
  nearest sample. A larger interval gives a smaller index and longer walks, the
  same trade as the sampled suffix array in an FM-index.
- An **anchor** is a reference node that most haplotypes in the region visit,
  one every `--anchor-spacing` bp along each reference path. The haplotype index
  contains a sample for every visit to an anchor node, so the samples at that
  node list every haplotype passing it, with the position and coordinate of each
  visit.
- A **bin** is `--stray-bin` bp of a reference path, 16,384 by default. For each
  bin, the haplotype index lists the nodes that a query's `context` reaches from
  the reference nodes of the bin, for a `context` up to `--stray-context`. The
  subgraph of a window consists of nodes listed for the bins the window touches,
  plus the snarls the query fills.
- A **stray row** records a stretch of one path that visits the nodes of a bin
  where a walk from one anchor to the next can miss it. A contig that visits no
  anchor has stray rows, and so have the stretch from a contig's end to its
  first anchor visit and a pass over nodes that the graph shares with a copy of
  the region far along the path. Each row gives the bin, the path, the stretch
  and the GBWT position of its first visit.

![graph.gbz.db lists the walks at each node by rank and names each path at the start of its walk. The haplotype index adds a sample every --interval bp that maps a position to a path, so a query names a walk from the next sample along it](img/haplotype-samples.svg)

A query that uses the `keep` option reads the anchor rows around the window, and
the node lists and stray rows of the bins the window touches.
`gbz-haplotype-index` starts the interval count of each path again at every
anchor visit, so the samples of all haplotypes fall near the same reference
positions, and a window shorter than `--interval` often contains none. The
anchor rows give the position and coordinate of each chosen haplotype at each
anchor, and the query walks each chosen haplotype from one anchor to the next,
then along its stray rows ([below](#keep)).

Table `HaplotypeSamples` contains the samples, `HaplotypeAnchors` lists the
anchor nodes, `HaplotypeBinNodes` lists the nodes of each bin, `HaplotypeStrays`
contains the stray rows, and `HaplotypeLengths` lists the length of each path.

The word "sample" also means an individual, such as `HG002` in a path name. The
rest of this page uses it for the index entry.

## Building the haplotype index

```bash
cargo install gbz-haplotype-index

gbz-haplotype-index --interval 16384 graph.gbz graph.gbz.db graph.haplotype-index.db

# without the GBZ, reading paths from the database
gbz-haplotype-index --interval 16384 --from-db graph.gbz.db graph.haplotype-index.db
```

Give `graph.gbz.db` with the GBZ. The tool checks that the two match, and reads
the top-level snarls from the database, which the stray rows need for a query
that fills snarls. `--overwrite` replaces an existing index. `--anchor-spacing`
defaults to 32,768 bp, and `--anchor-sample` places the anchors, bins and stray
rows on the paths of one reference sample. A query on a reference path without
them identifies every walk. `--anchor-spacing 0` writes none. The source is in
`tools/haplotype-index/`.

`--stray-context` sets the largest `context` that the node lists and stray rows
cover, 1,000 bp by default. A query with a larger `context` identifies every
walk, and `--stray-context 0` writes neither table. `--stray-bin` sets the bin
size, `--stray-bound` how far a walk runs past an anchor visit (32,768 bp), and
`--stray-gap` the distance between two visits that starts a new stray row (1,024
bp). `--reference-interval` sets the interval between samples along the
reference paths, and defaults to `--interval`.

Keep the default of sampling both orientations. `open` rejects an index built
with `--forward-only`, because about half the contigs in a graph like HPRC's run
reversed relative to the reference, and identifying a walk on one of those needs
reverse-orientation samples.

For the 10 GB HPRC v2.1 GRCh38 database, the index with 131,072 bp anchors on
GRCh38 and CHM13 is 8.1 GB. It holds 178.5 million samples, 5.4 million stray
rows, and node lists of 3.2 MB for 363,000 bins; the stray rows and node lists
add 0.2 GB. Building it from the GBZ takes 33 minutes on 20 threads and peaks at
19.5 GB of memory.

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
visit: the walks the same query without `keep` returns for those haplotypes.
`--keep` takes a sample or `sample#haplotype` and can repeat. A query that uses
the `keep` option needs the haplotype index, because the library looks up the
haplotype of each walk in it. With `haplotypes: 'distinct'`, the query merges
the identical walks of the haplotypes you kept, and each `weight` counts kept
haplotypes.

## How a query identifies walks

The library identifies walks by one of two routes, and both return the same
walks with the same GBWT positions and coordinates. The tests compare the two
routes record for record, and check each record by walking back through the GBWT
to a sample of its path. [Measured on HPRC](#measured-on-hprc-v21) lists the
comparisons on a whole graph.

![The sampled route extracts every walk through the window and follows each one to a sample. The keep route reads the rows at the anchor nodes around the window, then walks the chosen haplotypes from one anchor to the next. A strip of chr6 marks one anchor per 131,072 bp and the two that a window at 33,000,000 lies between](img/naming-routes-layout.svg)

- **Sampled**: the query extracts every walk through the window, then follows
  each walk to a sample to identify it. The cost grows with the number of
  haplotypes in the window.
- **Keep**: the query builds the same subgraph and finds the walks of the chosen
  haplotypes in it from the anchor rows and the stray rows. The cost grows with
  the number of haplotypes chosen.

A query that uses the `keep` option takes the keep route when the haplotype
index has node lists and stray rows, which `gbz-haplotype-index` 0.2.0 writes.
The keep route returns its walks when its checks show that they are complete
([below](#keep)). When a check fails, and with a haplotype index that an earlier
`gbz-haplotype-index` built, the library identifies every walk on the same
subgraph, then drops the haplotypes the predicate rejects.

![A query that uses the keep option takes the keep route; when a check of the keep route fails, and for every other query, the library identifies every walk](img/naming-routes.svg)

The flowchart source is [naming-routes.dot](img/naming-routes.dot); the
schematic is hand-written SVG.

### Sampled

1. The query extracts every walk through the nodes in the window.
2. `identifyPaths()` identifies each walk from a sample the walk passes inside
   the window. A walk with no sample there is followed past the window, for
   about four sampling intervals, and printed as `unknown#N` if that finds none.
3. For a query that uses the `keep` option, the library then drops the
   haplotypes the predicate rejects. It throws when a walk has no name, because
   that walk could belong to a haplotype you kept.

### Keep

1. The query builds the subgraph as the sampled route does, from the reference
   walk through the window, `context` bp around it and the snarls that `snarls`
   selects. A piece of a walk is a run of its positions whose nodes all lie in
   this subgraph.
2. The query reads the node lists of the bins the window touches, and checks
   that each node of the subgraph appears in one. A node that the query added by
   filling a snarl passes the check when the stray rows cover that snarl.
3. The query reads the anchors within 32 kb of those bins, plus one more on each
   side, and the rows at those anchors.
4. The query walks each chosen haplotype's stretch between visits to two
   adjacent anchors whose span meets the bins, from 32 kb before the first visit
   to 32 kb past the second, and records every piece on the way. When the span
   ends up to 32 kb outside the bins, the query walks 32 kb to each side of the
   visit at that end.
5. The query reads the stray rows of the bins, and walks each row of a chosen
   haplotype from its first visit to its last.
6. The query reads the samples on the nodes of the subgraph and checks that each
   sample of a chosen haplotype lies in a piece it found.
7. `extractPaths` keeps each piece in the orientation whose end nodes are
   canonical, and the walks of steps 4 and 5 find each piece in the forward
   orientation of its path. A reverse-orientation sample from step 6 gives a
   position in the other orientation of its piece. For a piece without one, the
   query walks on past the piece and reads the reverse-orientation samples at
   each node until one matches the coordinate, then walks from that sample into
   the piece.

`gbz-haplotype-index` writes the stray rows by the rule that steps 3 and 4 walk
by. For each bin it finds the nodes with the search that a query's `context`
runs, started from every reference node of the bin. It then follows every path
and tests each visit to a node against each bin that lists the node. Steps 3 and
4 reach a visit that lies in a stretch between visits to adjacent anchors whose
span meets the bin, and a visit within 32 kb along the path of an anchor visit
that lies up to 32 kb outside the bin. Every other visit goes into a stray row
of that bin. Step 2 checks that the subgraph consists of nodes of the bins the
window touches, so steps 4 and 5 pass every visit of a chosen haplotype to those
nodes.

A snarl that the query fills can add nodes that no bin lists. A path that passes
one of those nodes enters the snarl through one of its two boundary nodes, which
the bins list, so a walk passes that visit and the piece runs on into the snarl.
A path that lies inside the snarl from end to end has a stray row that names the
snarl, and the query walks that row when it filled the snarl.

The query identifies every walk when one of these holds, and
`gbz-base-query --stats` prints which:

- its `context` exceeds `--stray-context`;
- step 2 finds a node that no bin it touches lists;
- `snarls` is `'overlapping'` and the query filled a snarl;
- the query filled a snarl of more than 1,048,576 nodes, or any snarl when the
  index was built without `graph.gbz.db`;
- more than 32 chosen paths pass the anchors, because the sampled route took
  less time than the walks for 42 haplotypes;
- a walk goes farther outside the subgraph than eight anchor spacings plus twice
  the window's length;
- a sample in step 6 lies outside every piece;
- a twin in step 7 stays out of reach.

### Measured on HPRC v2.1

We ran the index described under
[Building the haplotype index](#building-the-haplotype-index) from local copies
of both files, with the scripts in `tools/validate/`.

Each query ran on the keep route and on the sampled route, and we compared the
GFA and the alignment records of the two:

| Windows                                                                                                                                                          | Queries | Keep route | Identified every walk | Differ |
| ---------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------: | ---------: | --------------------: | -----: |
| 1,500 of 300 bp to 100 kb: random, segmental duplications, chromosome ends, acrocentric arms, CHM13, unplaced contigs, anchors that a path visits twice          |   9,000 |      8,986 |                    14 |      0 |
| 1,330 aimed at the stray rows: the bin of a stray row, a path inside one snarl, the ends of an edge over 16 kb of reference, bin boundaries, anchors, 150-500 kb |  15,780 |     15,772 |                     8 |      0 |
| 18 at KIV-2, AMY1, MHC, SMN, GSTM1, C4, CFH, and seven places where an earlier route dropped pieces                                                              |     216 |        216 |                     0 |      0 |
| the first and third rows again, keeping every haplotype                                                                                                          |   3,072 |      3,070 |                     2 |      0 |

The windows of the first row ran with `context` 100 and no snarls, and with
`context` 1000 and contained snarls. Those of the second and third also ran with
`context` 0 and contained snarls, and with `context` 1000 and no snarls. Each
kept one haplotype, one sample, and eight haplotypes. The last row kept every
haplotype, with the limit of 32 chosen paths lifted, so it compares every walk
in each window. The 24 queries that identified every walk had more than 32
chosen paths at the anchors (18) or a walk past its cap (6). The route of 4.1.0,
which placed each haplotype from its anchor visits and the samples in the
window, differed in 152 of the first row's queries and dropped 585 pieces.

`gbz-truth` lists the pieces that every path leaves in a subgraph, from the GBZ
alone. The sampled route returned those pieces, with the same path and
coordinates, in each of the 8,260 subgraphs of the first two rows: 21,653,445
pieces.

`test/fuzz/` generates graphs with inversions, duplications, contigs split into
fragments and contigs that visit no anchor, and reads the pieces from the GFA.
Over 18,500 graphs, 582,634 sampled queries and 7,322,972 queries that use the
`keep` option returned those pieces, 6,974,349 of the second on the keep route.

An audit then set the index options against the route, with samples at anchors
and path ends only, bins of 1 bp, bounds under the spacing, and windows aimed at
the bound around each anchor. The keep route answered 1,706,363 queries in seven
such runs, and none differed from the GFA. The audit found one omission by
reasoning, at an anchor whose node starts exactly `bound` past the window's
bins, and fixed it with a fixture (`test/data/anchor-at-bound.gfa`). The trigger
needs a node longer than half the anchor spacing, and the longest node in HPRC
v2.1 is 342 bp.

With both files on local disk and their pages in memory, the median query took
102 ms on the keep route and 169 ms on the sampled route over the random
windows, 117 ms and 289 ms in segmental duplications, and 144 ms and 1,404 ms at
unplaced contigs. On a window under 3 kb beside an anchor or a bin boundary, the
keep route took 62-93 ms and the sampled route 19-25 ms, because the keep route
walks each chosen haplotype from one anchor to the next, 131 kb apart in this
index, whatever the window's length.
