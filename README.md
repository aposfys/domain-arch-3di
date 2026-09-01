# domain-arch-3di
Does Foldseek 3Di structural encoding change inferred protein domain rearrangement events, or only tree topology?

[![CI](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

58 human proteins across two clades — 29 MFS transporters (IPR020846), where domain rearrangement is claimed to be most active, and 29 trypsin-domain serine proteases (IPR001254) as a globular control. AlphaFold models, Foldseek 3Di, pLDDT masking, and the same alignment-free distance and neighbour-joining step applied to both alphabets.

### Both change

**Topology:** Robinson-Foulds 80, normalised **0.727**. About three-quarters of the bipartitions in one tree are absent from the other.

**Events:** read off sister pairs, which is where two extant architectures can be compared without reconstructing an ancestor.

| | Sequence tree | Structure tree |
| --- | ---: | ---: |
| Cherries | 21 | 18 |
| no change | 10 | 13 |
| complex | 7 | 3 |
| terminal addition | 3 | 1 |
| internal deletion | 1 | 1 |
| **Rearrangement events** | **11** | **5** |

Only 9 sister pairs are shared between the two trees (Jaccard 0.30). **So it is not the case that 3Di reshuffles the phylogeny while leaving the evolutionary story intact** — the inferred rearrangements move too, and here the structural characters imply less than half as many.

Of the 9 cherries both trees found, 0 were classified differently. That has to be zero and is reported rather than assumed: a non-zero value would mean the event classifier depends on something beyond the architectures themselves.

### Running it

```
make install
export FOLDSEEK_BIN=/path/to/foldseek
python3 -m domarch.cli fetch --per-clade 30       # UniProt + InterPro
python3 -m domarch.cli analysis --per-clade 30    # AlphaFold, 3Di, both trees
make test
```

### Two things that would have gone wrong quietly

- **AlphaFold model versions move.** The `v4` URL template that was current when this repo was designed now returns `NoSuchKey`, and a hardcoded template fills the cache with 127-byte XML error documents that Foldseek then fails on for an unrelated-looking reason. Download URLs are resolved through the API, and a payload under 1 kB is rejected as not-a-model.
- **A masked residue is not a letter.** pLDDT masking replaces low-confidence residues with `X`; any *k*-mer window touching one is dropped rather than counted, because counting `X`-containing k-mers turns disordered regions into their own signal — the exact artefact the masking exists to remove. A test pins it.

Proteins whose AlphaFold model length disagrees with their UniProt sequence are excluded rather than aligned by position, and proteins below 50% confident residues are excluded as a first-class step with the exclusions reported.

### Limitations

Both alphabets go through an identical 3-mer cosine distance, which removes the confound of comparing two substitution matrices that were never calibrated against each other — at the cost of resolution, so these trees are weaker than a model-based inference and the RF distance is noisier. Around twenty cherries per tree is enough to show the event sets differ and not enough to estimate rates. These are paralogues within one species, so the result is about how characters change an inference, not about the evolution of these families.

### Layout

```
src/domarch/
  data.py           UniProt clades and InterPro domain architectures
  structure.py      AlphaFold retrieval, 3Di encoding, pLDDT masking
  architecture.py   rearrangement event classification
  trees.py          k-mer distances, neighbour joining, splits, cherries
  compare.py        event-level comparison between the two trees
  analysis.py       the whole run
  report.py         results rendering
  cli.py            fetch / analysis / report
```

27 tests, none needing a network or a structure.

### More

- [Analysis: what was done, and why it was done that way](ANALYSIS.md)
- [Full results](results/RESULTS.md)
- [Where the field disagrees with itself, and the traps this pipeline avoids](docs/DESIGN.md)
