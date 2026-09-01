# domain-arch-3di
Does Foldseek 3Di encoding change inferred domain rearrangement events, or only tree topology?

[![CI](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

58 human proteins across two clades — 29 MFS transporters (IPR020846), where
rearrangement is claimed to be most active, and 29 trypsin-domain serine
proteases (IPR001254) as a globular control. AlphaFold models, Foldseek 3Di,
pLDDT masking, and the same alignment-free distance and neighbour-joining step
applied to both alphabets.

```
make install
export FOLDSEEK_BIN=/path/to/foldseek
python3 -m domarch.cli fetch --per-clade 30       # UniProt + InterPro
python3 -m domarch.cli analysis --per-clade 30    # AlphaFold, 3Di, both trees
make test                                         # 27 tests, no network, no structure
```

### Both change

**Topology:** Robinson-Foulds 80, normalised **0.727** — about three-quarters of
the bipartitions in one tree are absent from the other.

**Events**, read off sister pairs, where two extant architectures can be compared
without reconstructing an ancestor:

| | Sequence tree | Structure tree |
| --- | ---: | ---: |
| Cherries | 21 | 18 |
| no change | 10 | 13 |
| complex | 7 | 3 |
| terminal addition | 3 | 1 |
| internal deletion | 1 | 1 |
| **Rearrangement events** | **11** | **5** |

Only 9 sister pairs are shared between the trees (Jaccard 0.30). **So 3Di does not
merely reshuffle the phylogeny while leaving the evolutionary story intact** — the
inferred rearrangements move too, and the structural characters imply less than
half as many.

Of the 9 cherries both trees found, 0 were classified differently. That has to be
zero and is reported rather than assumed: a non-zero value would mean the event
classifier depends on something beyond the architectures themselves.

### Limitations

Both alphabets go through an identical 3-mer cosine distance, which removes the
confound of comparing two uncalibrated substitution matrices — at the cost of
resolution, so these trees are weaker than a model-based inference and the RF
distance is noisier. Around twenty cherries per tree is enough to show the event
sets differ and not enough to estimate rates. These are paralogues within one
species, so the result is about how characters change an inference, not about the
evolution of these families.

### More

- [Analysis](ANALYSIS.md) — what was done and why, including two failure modes that would have gone wrong quietly
- [Results](results/RESULTS.md) — full results
- [Design](docs/DESIGN.md) — where the field disagrees with itself, the layout, and the traps this avoids
