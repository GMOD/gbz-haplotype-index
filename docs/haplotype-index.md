# Naming haplotypes

Upstream gbz-base cannot tell which haplotype a subgraph walk belongs to, so it
prints `unknown#N`. This package names each walk using a haplotype index: side
tables that a small Rust tool, built on the unmodified upstream crates, writes
into an existing database or into a separate companion database.

## Building the index

```
cd tools/haplotype-index && cargo build --release
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 graph.gbz graph.gbz.db
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 --from-db graph.gbz.db
```

The second form walks the paths through the database's own node records, so it
works on a database whose GBZ is no longer at hand. Both forms write identical
tables. Walking a GBZ uses every core, or the number `--threads` sets.

## A companion database, for graphs you do not host

With `--output index.db` the tool writes the same tables into a standalone
companion database, and the reader opens the two side by side:

```
./target/release/gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --output graph.haplotype-index.db graph.gbz
gbz-base-query https://host/graph.gbz.db --haplotype-index https://host/graph.haplotype-index.db ...
```

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
```

A companion adds haplotype names to a database someone else publishes, without
rehosting it. HPRC publishes `hprc-v2.1-mc-grch38.gbz.db` (10 GB) beside its
graphs, and the companion for it is built from the 5 GB GBZ. The companion
records the graph's path and node counts, and the reader refuses a companion
built for a different graph.

## Tables

`HaplotypeSamples` holds one GBWT position every `--interval` bp along every
path, in both orientations, with the path handle and the forward coordinate of
that node. `HaplotypeLengths` holds each path's length.

`HaplotypeAnchors` names one anchor node per multiple of `--anchor-spacing`
along every path the database indexes for random access, which are its reference
paths. For the multiple 0 the anchor is the path's first node. For every other
multiple it is the node with the most GBWT positions among those overlapping the
half spacing before the multiple. Most haplotypes of that stretch pass through
such a node, whereas the node containing the multiple itself can be a rare
allele in a variable region. The tool also writes every path's visit through an
anchor node to `HaplotypeSamples`, in both orientations, so the rows at one
anchor node list every haplotype passing that point of the reference, each with
its own coordinate.

On the HPRC v2.1 graph at the default 131,072 bp spacing, the tool writes 45,557
anchors over the 292 GRCh38 and CHM13 paths, and the companion grows from 7.0 to
about 7.9 GB. `--anchor-spacing 0` writes no anchors. The `Tags` table records
the spacing and the rule, and the upstream `query` binary keeps working on an
augmented database.

The index needs both orientations because about half the contigs of a graph like
HPRC's are stored against their reference, and a forward-only index holds no
sample for a walk on one of those. `GBZBase.open` throws `ForwardOnlyIndexError`
for an index written with `--forward-only`, which its
`haplotype_index_orientations` tag records as `forward`, instead of returning
half-named results.

## The sampled walk

A query without `keep` names its walks with `subgraph.identifyPaths()`. The
method loads the samples for the window's nodes with one index scan per run of
consecutive node ids; a window whose nodes sit in far-apart id ranges, as a
tandem repeat's do, takes one scan per range. It then chains each haplotype's
fragments to the next through the private nodes between them. When a chain
contains no sample inside the window, `identifyPaths()` follows its path onward
until it reaches one, up to four intervals past the last fragment it linked.

`identifyPaths()` fills in the named half of each record: the PanSN name, the
interval in that contig's coordinates, and the path handle. The range queries
call it when the database has the tables. The lower-level queries leave it to
the caller, before `alignments()` or `toSubgraphJson({ names: 'resolved' })`. On
the command line, `--resolve` and `--alignments` run it.

## The anchored walk

A query with `keep`, on a companion that carries anchors, takes the anchored
route, which `--stats` reports by name. The reader:

1. looks up the anchor for the multiple of the spacing at or before the window;
2. walks the reference from that node to a little past the window, to learn
   which nodes are the reference's and at what offsets;
3. reads the rows at the anchor node;
4. walks each row whose path the predicate accepts forward with `lf()`, from the
   row's position through the window.

The reader takes each walk's identity and coordinate from its anchor row, so the
route needs no chaining and no index scan, and nothing outside the kept set is
extracted or named.

A path through a duplicated stretch has several rows at the anchor node, one per
visit, and only one visit leads into the window. The reader tries first the
visit whose row sits nearest the reference's row, since GBWT orders rows by the
sequence before them, and stops at the first that reaches the window. The reader
therefore skips HG01109#1's second amylase copy, where walking every visit to
the bound would take 30,000 steps.

A kept contig with no row at the anchor, because it bypasses that node or starts
inside the window, is found as in the sampled route, from its samples on the
window's nodes. The reader follows one sample in the orientation that runs with
the reference back to the reference before the window, and starts the walk
there. A walk in the other orientation would leave the window backwards.

When a walk cannot be completed, the reader falls back to the sampled route for
the whole window, and `--stats` reports why. On a companion without
`HaplotypeAnchors`, and for a query with `haplotypes: 'distinct'`, the reader
takes the sampled route and then trims to the kept set.
[performance.md](performance.md#keeping-a-set-of-haplotypes) measures both
routes.

## Cutting to a set of haplotypes

`keepHaplotypes(name => ...)` narrows an identified subgraph to the reference
walk and the walks whose `PathName` the predicate accepts, and drops every node
only the other walks visited. The range queries' `keep` option calls it, and on
the command line it is `--keep SAMPLE` or `--keep SAMPLE#HAP`, repeatable. The
result is the cut for that set of haplotypes, ready to draw as it is.

A named walk in GFA or JSON output lists its steps in the haplotype's own
direction, so its `start..end` and its steps agree, as the GFA W line requires.

## Verification

The tests check every named fragment against an independent backward walk
through the bidirectional GBWT to the path's recorded start position.
