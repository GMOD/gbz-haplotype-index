# Naming haplotypes

Upstream gbz-base cannot say which haplotype a subgraph path belongs to, so it
emits `unknown#N`. This package adds that with two side tables that a small Rust
tool writes into an existing database, built on the unmodified upstream crates.

## Building the index

```
cd tools/haplotype-index && cargo build --release
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 graph.gbz graph.gbz.db
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 --from-db graph.gbz.db
```

The second form walks the paths through the database's own node records, so a
database whose GBZ is no longer at hand can still be augmented; the two forms
write identical tables. Walking a GBZ uses every core (`--threads`).

## A companion database, for graphs you do not host

With `--output index.db` the tool writes the same tables into a standalone
companion database instead, and the reader opens the two side by side:

```
./target/release/gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --output graph.haplotype-index.db graph.gbz
gbz-base-query https://host/graph.gbz.db --haplotype-index https://host/graph.haplotype-index.db ...
```

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
```

This is how a database someone else publishes gets haplotype names without
anyone rehosting it: HPRC publishes `hprc-v2.1-mc-grch38.gbz.db` (10 GB) beside
its graphs, and the companion for it is built from the 5 GB GBZ. The companion
records the graph's path and node counts and the reader refuses one built for a
different graph.

## What the tables hold

`HaplotypeSamples` holds one GBWT position every `--interval` bp along every
path in both orientations, with the path handle and the forward coordinate of
that node, and `HaplotypeLengths` holds each path's length.

`HaplotypeAnchors` names one anchor node per multiple of `--anchor-spacing`
along every path the database indexes for random access (its reference paths):
the path's first node for the multiple 0, and otherwise the node with the most
GBWT positions among those overlapping the half spacing before the multiple, so
it is a node most haplotypes of that stretch pass rather than whichever node
happens to contain the multiple, which in a variable region can be a rare
allele. Every path's visit through an anchor node is written to
`HaplotypeSamples` as well, in both orientations, so the rows at one node list
every haplotype passing that point of the reference with its own coordinate.

On the HPRC v2.1 graph at the default 131,072 bp that is 45,557 anchors over the
292 GRCh38 and CHM13 paths, and the companion grows from 7.0 to about 7.9 GB.
`--anchor-spacing 0` writes none, and the `Tags` table records the spacing and
the rule. The upstream `query` binary keeps working on the augmented database.

Both orientations are needed: about half the contigs of a graph like HPRC's are
stored against their reference, and a walk of one of those meets no sample from
a forward-only index. `GBZBase.open` refuses an index the tool wrote with
`--forward-only` (its `haplotype_index_orientations` tag says `forward`) with
`ForwardOnlyIndexError`, so the half-named result never reaches a caller.

## The sampled walk

At query time `subgraph.identifyPaths()` loads the samples for the window's
nodes with one index scan per run of consecutive node ids (a window whose nodes
sit in far-apart id ranges, as a tandem repeat's do, is not one scan across the
gap), chains each haplotype's fragments to the next through the private nodes
between them, and for a chain that met no sample inside the window walks on
until it finds one, up to four intervals past the last fragment it linked. This
is what fills in the `resolved` half of a record: PanSN name, haplotype interval
in that contig's coordinates, and the path handle.

`getAlignmentsForRange` and `getSubgraphForRange` run it for you when the
database has the tables; on the lower-level path you call it yourself before
`alignments()` or `toSubgraphJson({ names: 'resolved' })`. On the command line,
`--resolve` and `--alignments`.

## The anchored walk

A query with `keep` on a companion that carries anchors takes a different route,
and `--stats` reports it as the anchored walk. The reader looks up the anchor
for the multiple of the spacing at or before the window, walks the reference
from that node to a little past the window to learn which nodes are the
reference's and where, reads the rows at the anchor node, and for each row whose
path the predicate wants walks that path with `lf()` from its own position
through the window, so its identity and coordinate come from the row and no
chain walk or index scan is needed.

A path through a duplicated stretch has several rows at the anchor node, one per
visit, and only one visit is followed by the window: the reader tries the visit
whose row sits nearest the reference's own first (GBWT rows are ordered by the
sequence before them) and stops at the first that goes through, so HG01109#1's
second amylase copy costs nothing where it used to walk 30,000 steps to the
bound.

A wanted contig with no row at the anchor, because it bypasses that node or
starts inside the window, is found the way the sampled route finds every
haplotype, from its per-path samples on the window's nodes: one sample of the
orientation that runs with the reference (both are indexed, and the other would
walk out of the window backwards) is followed back to the reference before the
window and the walk starts there.

When a walk cannot be completed, the whole window falls back to the sampled
route and `--stats` says why. A companion without `HaplotypeAnchors` takes the
sampled route and trims it, as does `haplotypes: 'distinct'`. What each route
costs is measured in
[performance.md](performance.md#keeping-a-set-of-haplotypes).

## Cutting to a set of haplotypes

`keepHaplotypes(name => ...)` narrows an identified subgraph to the reference
walk and the walks whose PanSN name the predicate accepts, dropping every node
only the other walks visited; it is what the range queries' `keep` option calls,
and on the command line `--keep SAMPLE` or `--keep SAMPLE#HAP` (repeatable).
This is how a cut for a chosen set of haplotypes is written once and drawn as it
is.

A named walk in GFA or JSON output lists its steps in the haplotype's own
direction, whichever twin of the walk the extraction kept, so `start..end` and
the steps agree as the W line spec requires.

## Verification

The tests check every resolved fragment against an independent backward walk
through the bidirectional GBWT to the path's recorded start position.
