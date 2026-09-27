# Naming haplotypes

Upstream gbz-base prints each walk in a window as `unknown#N`. This package
names the walks using a haplotype index: extra tables that the Rust tool in
`tools/haplotype-index/` writes into the database, or into a separate companion
file.

## Building the index

```
cd tools/haplotype-index && cargo build --release

# into the database
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 graph.gbz graph.gbz.db

# into the database, reading paths from the database itself
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 --from-db graph.gbz.db

# into a companion file
./target/release/gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --output graph.haplotype-index.db graph.gbz
```

A companion adds names to a database someone else hosts. Open the two together:

```ts
const db = await GBZBase.open(new RemoteFile(graphUrl), {
  haplotypeIndex: new RemoteFile(indexUrl),
})
```

```
gbz-base-query https://host/graph.gbz.db --haplotype-index https://host/graph.haplotype-index.db ...
```

For the 10 GB HPRC v2.1 GRCh38 database, the companion is 7.9 GB. It records the
graph's path and node counts, and `open` checks them against the graph.

## Tables

- `HaplotypeSamples`: a GBWT position every `--interval` bp along every path, in
  both orientations, with the path and its coordinate there.
- `HaplotypeLengths`: the length of each path.
- `HaplotypeAnchors`: one node every `--anchor-spacing` bp along each reference
  path, chosen as the node most haplotypes pass through near that offset. Every
  haplotype's visit to an anchor node is also a sample, so an anchor lists every
  haplotype passing that point.

`open` rejects an index built with `--forward-only` with
`ForwardOnlyIndexError`. About half the contigs in a graph like HPRC's are
stored reversed relative to the reference, and walks on those need
reverse-orientation samples.

## Two ways to name walks

**Sampled.** `identifyPaths()` reads the samples on the window's nodes, links
each haplotype's fragments into chains, and names each chain from a sample. A
chain with no sample in the window is followed past it, up to four sampling
intervals. The range queries run this for you; the lower-level queries leave it
to you, and the CLI runs it with `--resolve` or `--alignments`.

**Anchored.** A query with `keep` on a companion with anchors reads the samples
at the anchor before the window and walks only the chosen haplotypes forward
from there through the window. The cost scales with the number of chosen
haplotypes. If a walk fails, the query falls back to the sampled route, and
`--stats` reports the reason.
[performance.md](performance.md#keeping-a-set-of-haplotypes) compares the two.

## Keeping a set of haplotypes

`keepHaplotypes(predicate)`, which the `keep` option and `--keep SAMPLE[#HAP]`
call, reduces a named subgraph to the reference, the accepted walks and the
nodes they visit.

The tests check every named fragment by walking backward through the GBWT to the
path's recorded start.
