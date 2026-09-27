# Naming haplotypes

Upstream gbz-base names only the query path and prints every other walk as
`unknown#N`. This package names the walks using a haplotype index, a set of
extra tables that the Rust program `gbz-haplotype-index` writes into the
database or into a separate companion file. Build the index once per graph, then
pass it to the library or the command line.

## Building the index

The source is in `tools/haplotype-index/` and needs a Rust toolchain.

```bash
cd tools/haplotype-index && cargo build --release

# into the database
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 graph.gbz graph.gbz.db

# into the database, reading paths from the database itself
./target/release/gbz-haplotype-index --interval 4096 --anchor-spacing 131072 --from-db graph.gbz.db

# into a companion file
./target/release/gbz-haplotype-index --interval 16384 --anchor-spacing 131072 --output graph.haplotype-index.db graph.gbz
```

A companion file names the walks of a database someone else hosts. In the
library, open the two together:

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

For the 10 GB HPRC v2.1 GRCh38 database, the companion is 7.9 GB. The companion
records the graph's path and node counts, and `open` checks them against the
graph.

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

## Two ways to name walks

**Sampled.** `identifyPaths()` reads the samples on the window's nodes, links
each haplotype's fragments into chains, and names each chain from a sample. When
a chain has no sample in the window, `identifyPaths()` follows it past the
window, up to four sampling intervals. The range queries call `identifyPaths()`
for you, and the lower-level queries leave the call to you. On the command line,
`--resolve` and `--alignments` run it.

**Anchored.** A query with `keep` on a companion with anchors reads the samples
at the anchor before the window and walks only the chosen haplotypes forward
from there through the window. The cost scales with the number of chosen
haplotypes. If a walk fails, the query falls back to the sampled route, and
`gbz-base-query --stats` reports the reason.
[performance.md](performance.md#keeping-a-set-of-haplotypes) compares the two.

## Keeping a set of haplotypes

The `keep` option and `--keep SAMPLE[#HAP]` call `keepHaplotypes(predicate)`,
which reduces a named subgraph to the reference, the accepted walks and the
nodes they visit.

The tests check every named fragment by walking backward through the GBWT to the
path's recorded start.
