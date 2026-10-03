# The haplotype index

Upstream gbz-base prints each walk other than the query path as `unknown#N`.
[@gmod/gbz-base](https://github.com/GMOD/gbz-base-js) reports the sample,
haplotype and contig of every walk, looked up in a haplotype index: a sidecar
SQLite file that this repository's program,
[`gbz-haplotype-index`](https://crates.io/crates/gbz-haplotype-index)
(`cargo install gbz-haplotype-index`), writes beside the graph database. The
graph database stays as `gbz-base construct` wrote it, so a haplotype index also
works with a database someone else hosts.

Building the haplotype index happens once per graph and covers every haplotype.
Identifying walks happens on every query that reads the index. A caller who
wants a subset of the haplotypes passes the `keep` option with the query
([below](#querying-a-subset-of-the-haplotypes)), with no rebuild.

## Samples and anchors

A path is one contig of one haplotype, named `sample#haplotype#contig`. The GBWT
stores each path as a sequence of oriented nodes, once in each orientation. A
position is a node record plus the rank of one visit among all visits to that
node, so each position belongs to exactly one path in one orientation, and each
node record maps every position to the next one along its path. A GBWT file also
carries document array samples that map a position back to its path, but gbwt-rs
does not read them, so the database names a path only at its two start positions
in `Paths`, which can be a whole chromosome away from a query window. The
haplotype index is that sampling, kept beside the database. It records the path
at positions along the way:

- A **sample** records a position, the path through it, the orientation and the
  coordinate along that path. `gbz-haplotype-index` writes one every
  `--interval` bp along each path, in both orientations, plus one at the start
  and end of each path. To identify a walk, the library follows the walk to the
  nearest sample. A larger interval gives a smaller index and longer walks, the
  same trade as the sampled suffix array in an FM-index.
- An **anchor** is a reference node that most haplotypes in the region visit,
  one every `--anchor-spacing` bp along each reference path. A reference path
  here is one gbz-base indexes for random access, which are the paths of the
  samples in the GBWT tag `reference_samples` and the generic paths, whose
  sample is `_gbwt_ref`. The haplotype index contains a sample for every visit
  to an anchor node, so the samples at that node list every haplotype passing
  it, with the position and coordinate of each visit.
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

![graph.gbz.db lists the visits at each node by rank and names each path at its start position. The haplotype index adds a sample every --interval bp that maps a position to a path, so a query names a walk from the next sample along it](img/haplotype-samples.svg)

`gbz-haplotype-index` starts the interval count of each path again at every
anchor visit, so the samples of all haplotypes fall near the same reference
positions, and a window shorter than `--interval` often contains none. A query
that uses the `keep` option therefore reads the anchor rows, node lists and
stray rows around the window instead ([below](#keep)).

Table `HaplotypeSamples` holds the samples in the b-tree of their key, the node
and the rank of the visit. `HaplotypeAnchors` names each anchor node and holds,
as one blob per anchor, the samples of every visit through it. `HaplotypeBins`
holds each bin's node list beside its stray rows, and `HaplotypeLengths` the
length of each path. The tables with blobs keep SQLite's rowid and a key index,
their rows in key order, because an index b-tree copies whole rows into its
interior pages and a blob would leave each page a few cells. The tag
`haplotype_index_format` names this layout, 3. The library also reads format 2,
which gbz-haplotype-index 0.2 wrote: the samples in a rowid table, the anchors
without their visits, and the node lists and stray rows in two tables. A
keep-route query read that layout in about twenty single-block requests, most
of them b-tree descents four levels deep; format 3 answers the same query in
two short range scans ([measured below](#requests-per-query)).

### The overview

A view of megabases, or of a whole chromosome, is too much graph to read: 3 Mb
of chr22 costs about 70 requests and 17 MB, and whole chr6 timed out after
4,200. For those views the index carries an overview, which the library reads
with `haplotypeOverview` and `gbz-base-query --overview` in about 7 requests
whatever the window. For each reference sample with anchors, the indexer walks
every path once more. A path's visits to the sample's reference nodes place it
on a reference path; the stretch between two consecutive visits is an
excursion, which replaces the reference bp between the two nodes with the
path's own bp. The path covers the reference from the first node to the second,
and an excursion of `--overview-sv` bp or more (50) in either measure marks the
bins it spans as variant. A contig reads the reference in either direction, so
two visits in the same orientation that step on in that direction are an
alignment; a turn, a step back, a jump to another reference path or one of
over 10 Mb marks both nodes' bins as variant instead.

The paths of one sample and phase make one haplotype, and per bin of
`--overview-bin` bp (4,096) it gets one of four classes: absent when none of
its paths covers the bin, variant when one marked it, partial when they cover
less than nine tenths of it, and reference-like otherwise. Coarser levels, a
factor of four apart, sum the coverage and the marks of their bins, and a
variant cell carries a bucket of the marks (1, 2 to 3, 4 to 15, 16 or more), so
haplotypes still differ at a level where nearly every bin holds some structural
difference. `HaplotypeOverviewClasses` holds the cells, four bits each, in rows
of `--overview-chunk` bins (256) per reference path and level;
`HaplotypeOverviewBins` counts the haplotypes of each class and the excursions
per bin; `HaplotypeOverviewRows` names the haplotypes. On chr22 the overview
of 464 haplotypes takes 15 s per reference sample and 6.8 MB. Whole chr22 at
16 kb bins is 3,102 bins and 1.4 MB in 11 index requests, 4 of them once per
session. The build walks every path once more per reference sample and holds
8 bytes per node plus 3 bytes per haplotype and bin, about 2 GB more on HPRC
v2.1 than the figures above.

The word "sample" also means an individual, such as `HG002` in a path name. The
rest of this page uses it for the index entry, the sense of the GBWT's own
document array samples.

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
[`src/`](../src).

`--page-size` sets the SQLite page size, 65,536 by default, so that one page is
the 64 KiB block a remote reader fetches and every table is two or three levels
deep. `--stray-context` sets the largest `context` that the node lists and
stray rows cover, 1,000 bp by default. A query with a larger `context` identifies every
walk, and `--stray-context 0` writes neither table. `--stray-bin` sets the bin
size, `--stray-bound` how far a walk runs past an anchor visit (32,768 bp), and
`--stray-gap` the distance between two visits that starts a new stray row (1,024
bp). `--reference-interval` sets the interval between samples along the
reference paths, and defaults to `--interval`.

Keep the default of sampling both orientations. `open` rejects an index built
with `--forward-only`, because about half the contigs in a graph like HPRC's run
reversed relative to the reference, and identifying a walk on one of those needs
reverse-orientation samples.

For the 10 GB HPRC v2.1 GRCh38 database, the format 2 index with 131,072 bp
anchors on GRCh38 and CHM13 is 8.1 GB. It holds 178.5 million samples, 5.4
million stray rows, and node lists of 3.2 MB for 363,000 bins; the stray rows
and node lists add 0.2 GB. Building it from the GBZ takes 30-33 minutes on 20
threads and peaks at 19.9 GB of memory. On the chr22 part of the same graph,
1,131 paths and 3.1 million nodes, the same options give a 98 MB index in
format 2 and a 58 MB one in format 3, which drops the rowid table's second copy
of every sample key.

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
nodes ([proof](#proof-that-the-keep-route-finds-every-piece)).

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

### Proof that the keep route finds every piece

When the keep route returns without falling back, it has recorded every piece
that each chosen path leaves in the subgraph. The proof below follows the code
on both sides,
[`src/chosenPaths.ts`](https://github.com/GMOD/gbz-base-js/blob/main/src/chosenPaths.ts)
in gbz-base-js and [`src/strays.rs`](../src/strays.rs), so a change to either
can be checked against it. The fuzzer tests the same claim, and the proof names
the comparisons it rests on. [proof.md](proof.md) states the same argument as
definitions, lemmas and a theorem, without the code identifiers.

**Notation.** R is the query's reference path, `spacing` the anchor spacing,
`bin` the bin length and `bound` the walks' reach (`--stray-bound`). The window
touches bins F to L of R, which span `[lo, hi)` with `lo = F × bin` and
`hi = (L + 1) × bin`. Bin b spans `[b_lo, b_hi)`, and `lo ≤ b_lo < b_hi ≤ hi`
for every b from F to L. S is the subgraph and I the nodes that snarl fills
added to it. P is a chosen path, and a visit is one step of P's forward walk,
placed at the offset along P where its node starts. `a(k)` is the offset along R
of the anchor of multiple k.

**Assumptions.** Each holds by construction in the indexer or the reader. The
library checks parts of assumptions 3 and 4 and none of the others:

1. The multiples with an anchor on R run from 0 without a gap, and
   `a(k) ≤ a(k + 1)`, with equality only when one node anchors both multiples
   (`mark_anchors` in `main.rs`). Every node is at least 1 bp long.
2. The index has a sample at every visit of every path to an anchor node (`walk`
   in `main.rs`), so the anchor rows list every such visit. Each forward sample,
   anchor and stray row records a visit's GBWT position and the offset along the
   path where its node starts, the quantity a walk's `left` tracks.
3. The indexer (`Snarls::find`) and the query (`fillBetween`) fill a snarl with
   the same search from the same chain link: from both boundary nodes inward,
   along edges on both sides of each node, never entering the two boundary
   nodes. In a GBZ the graph's edges are the GBWT's, so both searches follow the
   same edges, and one of those edges joins any two consecutive visits of a
   path. So the next visit after a visit to the filled region lies in the region
   or at a boundary node. The query uses snarl rows only when the index's tags
   say `modeled` and record the same number of chain links as the database, a
   count that does not compare the links themselves.
4. The index was built from this graph. When it opens the index, the library
   compares the path and node counts, which would miss a different graph with
   the same counts.

**Part 1: what a walk covers.** A walk covers a range of offsets along P when
every visit to S in that range lies in a piece the query recorded. A walk starts
at a GBWT position with that visit's offset (assumption 2) and adds each node's
length as it steps, so `left` stays the offset of the current visit.
`walkForward` steps on while a node starts at or before `until`, and
`walkBackward` steps back while the current node starts after `until`, so both
comparisons are inclusive. On a node of S each records the maximal piece through
that position and jumps to the piece's far end, and the visits it jumps over lie
in that piece. So each walk steps through a contiguous stretch of P, and for P's
anchor visits at offsets `u ≤ v`:

- a section walk covers `[u − bound, v + bound]`;
- a walk around the anchor visit at u covers `[u − bound, u + bound]`;
- the walk of a stray row covers `[path_start, path_end]` of the row.

A walk that meets an end of P stops early, and the range beyond holds no visits.
`markWalked` records the planned range widened to what the walk reached, so
every recorded range is covered, and the query skips a plan or a row only when
`covered` finds its range inside a recorded range of P.

**Part 2: the indexer's rule.** The indexer cuts P at its visits to the anchor
nodes of R's sample. A section is the stretch between two consecutive cuts that
are anchors of adjacent multiples k and k + 1 of R, where neither node anchors
another multiple of any reference path of R's sample. For a bin b of R, the
indexer counts a visit of the section as reached when one of these holds:

- (i) `a(k) < b_hi` and `a(k + 1) ≥ b_lo`;
- (ii) `a(k) ≥ b_hi`, `a(k) − b_hi ≤ bound`, and the visit lies within `bound`
  along P of the visit to the anchor at `a(k)`;
- (iii) `a(k + 1) < b_lo`, `b_lo − a(k + 1) ≤ bound`, and the visit lies within
  `bound` along P of the visit to the anchor at `a(k + 1)`.

The indexer tests a visit only against the sections it lies in, two when the
visit is itself a cut. A visit within `bound` of an anchor visit in the
neighbouring section therefore becomes a stray, which adds a row. The indexer
writes every visit to a node of b that none of these reaches into a stray row of
b, and a row runs from its first stray to its last.

![For a window of one bin b: in case (i) the anchors' span meets b and the query walks the section and bound past both visits; in case (ii) a(k) lies up to bound past b and the query walks bound to each side of its visit; case (iii) mirrors it before b; every other visit to a node of b is in a stray row. Below, the query reads every anchor from the last before lo − bound to the first past hi + bound, so an anchor at exactly hi + bound brings in the one after it](img/proof-reach.svg)

**Part 3: the query reads both anchors of every such section.** Each case of
part 2 gives `a(k) ≤ hi + bound` and `a(k + 1) ≥ lo − bound`. The query reads
every multiple from `lowest` to `highest`, moving `lowest` down while
`a(lowest) ≥ lo − bound` or `lowest` has no anchor, and `highest` up while
`a(highest) ≤ hi + bound`.

- If `k < lowest`, then `k + 1 ≤ lowest`, so
  `a(k + 1) ≤ a(lowest) < lo − bound`, a contradiction, or `lowest` is 0 and k
  cannot be below it.
- If `k + 1 > highest`, then `k ≥ highest`. Either `a(highest) > hi + bound`, so
  `a(k) > hi + bound`, a contradiction, or `highest` has no anchor, and by
  assumption 1 neither does `k + 1`.

So the query reads both anchors. The proof needs `≥` in the first loop and `≤`
in the second. With `<` in the second, an anchor at exactly `hi + bound` stopped
the loop before the multiple past it, which is the omission that `4c7f8f8`
fixed. The first loop never moves: `mark_anchors` chooses for multiple k a node
that starts before `k × spacing`, so `a(k) < k × spacing`, and the first read
already starts one multiple below `(lo − bound) / spacing`. The second loop
moves when an anchor's node starts far before its multiple. Mutating the first
loop's comparison therefore changes nothing, and no test can catch it.

The query's cuts of P are its visits to the anchors it read on R, a subset of
the indexer's cuts: R is one of the reference paths of its sample, so its anchor
nodes are cuts for the indexer, and the query reads the samples at both
orientations of each anchor node once and keeps the forward ones, which by
assumption 2 are all of P's visits there. Two cuts that are consecutive for the
indexer are therefore consecutive for the query, and a node that anchors one
multiple for the indexer anchors one among those the query read. So the query
sees every section of part 2 as a consecutive pair of anchor visits with
adjacent multiples.

**Part 4: every reached visit lies on a planned walk.** Take a section of part 2
with anchor visits at u ≤ v along P.

- Case (i): `a(k) < b_hi ≤ hi` and `a(k + 1) ≥ b_lo ≥ lo`, so the query plans
  the section walk, which covers the section by part 1.
- Case (ii): if `a(k) < hi`, then `a(k + 1) ≥ a(k) ≥ b_hi > lo` and the query
  plans the section walk, which covers `[u − bound, v + bound]` and with it
  every visit within `bound` of either anchor visit. Otherwise `a(k) ≥ hi` and
  `a(k) − hi ≤ a(k) − b_hi ≤ bound`, so the query plans a walk around the visit
  to the anchor at `a(k)`.
- Case (iii) is case (ii) reflected: the section walk when `a(k + 1) ≥ lo`, and
  otherwise a walk around the visit to the anchor at `a(k + 1)`.

**Part 5: every visit to S lies in a recorded piece.** Let p be a visit of P to
a node x of S.

- When x is not in I, the check of step 2 of the query found x in the node list
  of some bin b from F to L of R. The indexer tested p against b. If one of the
  cases of part 2 reached p, part 4 gives a planned walk that covers it.
  Otherwise p lies in a stray row of b for P, and the query reads the rows of
  bins F to L, keeps those of chosen paths and walks each.
- When x is in I, x lies in the region of a filled snarl with boundary nodes y
  and z. The query fills a contained snarl only when y and z were in S before
  any fill, so neither is in I, and the fill puts the whole region in S. By
  assumption 3, P leaves the region only through a visit to y or z. If P has
  such a visit next to its run through p, that visit lies in a recorded piece by
  the case above, and p lies in the same piece, since every node between them is
  in S. Otherwise all of P lies inside the snarl, and the indexer wrote a row
  naming the snarl in every bin that lists the lower boundary node. The query's
  check found that node in a bin from F to L, so the query reads the row and
  walks all of P.

![A filled snarl between boundary nodes y and z, both in the subgraph before any fill. A path that enters the region leaves it only through y or z, so its piece through the region holds a visit to y or z. A path inside the region from end to end has a row naming the snarl](img/proof-snarl.svg)

Every visit of P to S therefore lies in a recorded piece, and each recorded
piece is a maximal run of P in S, so the query has recorded every piece of P.

**What the proof leaves to the checks.** The walks find each piece in P's
forward orientation, and step 7 of the query finds the other orientation where
`extractPaths` keeps it, or falls back. `pieceThrough` extends a piece while the
next position is in S, so each recorded piece is a maximal run of P in S. The
proof uses neither the `context` fallback nor any agreement between the
indexer's search for a bin's nodes and the query's `context` search, because the
check of step 2 establishes what it needs. That agreement keeps the check from
failing. The sample check of step 6 never fails while the proof holds, and it
guards against an index built for another graph or a rule changed on one side
alone.

### Measured on HPRC v2.1

We ran the index described under
[Building the haplotype index](#building-the-haplotype-index) from local copies
of both files, with the scripts in
[`tools/validate/`](https://github.com/GMOD/gbz-base-js/blob/main/tools/validate).

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
`context` 0 and contained snarls, and with `context` 1000 and no snarls, except
the 30 windows of 150-500 kb, which ran the first two settings only. Each kept
one haplotype, one sample, and eight haplotypes. The last row kept every
haplotype, with the limit of 32 chosen paths lifted, so it compares every walk
in each window. The 24 queries that identified every walk had more than 32
chosen paths at the anchors (18) or a walk past its cap (6).

`gbz-truth` lists the pieces that every path leaves in a subgraph, from the GBZ
alone. The sampled route returned those pieces, with the same path and
coordinates, in each of the 8,260 subgraphs of the first two rows: 21,653,445
pieces.

[`test/fuzz/`](https://github.com/GMOD/gbz-base-js/blob/main/test/fuzz) in
gbz-base-js generates graphs with inversions, duplications, contigs split into
fragments and contigs that visit no anchor, and reads the pieces from the GFA.
Over 18,500 graphs, 582,634 sampled queries and 7,322,972 queries that use the
`keep` option returned those pieces, 6,974,349 of the second on the keep route.

An audit then set the index options against the route, with samples at anchors
and path ends only, bins of 1 bp, bounds under the spacing, and windows aimed at
the bound around each anchor. The keep route answered 1,706,363 queries in seven
such runs, and none differed from the GFA. The audit found one omission by
reasoning, at an anchor whose node starts exactly `bound` past the window's
bins, and fixed it with a fixture
([`test/data/anchor-at-bound.gfa`](https://github.com/GMOD/gbz-base-js/blob/main/test/data/anchor-at-bound.gfa)).
The trigger needs a node longer than half the anchor spacing, 65,536 bp in this
index, and the longest node in HPRC v2.1 is 1,024 bp.

The fuzzer now builds each database with top-level chains from a vg distance
index, so its queries fill snarls. Over 2,050 more graphs, 363,623 queries that
use the `keep` option returned the pieces from the GFA, 60,468 of them after
filling a snarl, and the two routes never differed.

After the last change to either route, at `aafb2bc`, we ran the comparisons in
the table and `gbz-truth` again, and every count came out the same.

With both files on local disk and their pages in memory, the median query took
102 ms on the keep route and 169 ms on the sampled route over the random
windows, 117 ms and 289 ms in segmental duplications, and 144 ms and 1,404 ms at
unplaced contigs. On a window under 3 kb beside an anchor or a bin boundary, the
keep route took 62-93 ms and the sampled route 19-25 ms, because the keep route
walks each chosen haplotype from one anchor to the next, 131 kb apart in this
index, whatever the window's length.

### Requests per query

Over HTTP, a query's cost is the number of range requests more than the bytes,
since the library issues the reads of one b-tree descent one after another.
[`tools/requests/`](https://github.com/GMOD/gbz-base-js/blob/main/tools/requests)
in gbz-base-js counts them. The table gives the median over three random chr22
windows per size, cold cache, served from local copies of the chr22 database
(196 MB) and its index, built with the HPRC options above. "Before" is
@gmod/gbz-base 6.0.1 reading format 2; "after" is the next version reading
format 3, with its parallel reads and read-ahead. The sampled route returns all
464 haplotypes; keep1 keeps HG002's two.

| window    | route   | before: requests (graph + index) | before: MB | after: requests (graph + index) | after: MB |
| --------: | ------- | -------------------------------: | ---------: | ------------------------------: | ---------: |
|       300 | sampled |                       14 (8 + 6) |        0.9 |                      13 (8 + 5) |        0.9 |
|       300 | keep1   |                     28 (10 + 18) |        2.7 |                     23 (10 + 13) |        2.4 |
|    10,000 | sampled |                       13 (8 + 5) |        0.9 |                      13 (8 + 5) |        0.9 |
|    10,000 | keep1   |                     28 (11 + 17) |        2.3 |                     23 (11 + 12) |        2.0 |
|   100,000 | sampled |                       17 (8 + 9) |        1.4 |                      13 (8 + 5) |        1.2 |
|   100,000 | keep1   |                     31 (11 + 20) |        2.7 |                     23 (11 + 12) |        2.3 |
|   300,000 | sampled |                      29 (11 + 18) |        2.9 |                     16 (11 + 5) |        2.6 |
|   300,000 | keep1   |                      45 (14 + 31) |        4.5 |                     27 (14 + 13) |        3.7 |
| 1,000,000 | sampled |                      56 (12 + 44) |        7.1 |                     17 (12 + 5) |        6.0 |
| 1,000,000 | keep1   |                      70 (15 + 55) |        8.6 |                     28 (15 + 13) |        7.3 |
| 3,000,000 | sampled |                   307 (181 + 126) |       20.1 |                     70 (63 + 7) |       16.6 |
| 3,000,000 | keep1   |                   324 (188 + 136) |       21.2 |                     81 (66 + 15) |       17.7 |

Of a keep1 query's 12 or 13 index requests, 2 open the file, 1 reads the path's
length, 2 or 3 each read the bins' rows and the anchors' rows (a key index,
then a run of the table), in parallel, and 3 read the samples on the window's
nodes for the check after the walks, which now run while the walks do. The graph side is unchanged below 1 Mb: the
reference walk's node records come in one or two prefetched reads, and the
rest are nodes the context expansion and the identification reach outside
them. Past 1 Mb the pager used to refuse a prefetch over half its cache and the
walk read one block per request; it now fetches the range in chunks as the
walk approaches them. The 3 Mb windows that still cost 60 to 180 graph
requests lie in a region where consecutive reference nodes sit a megabyte
apart in the node table, so each reference node is its own read.
