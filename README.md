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

### Prior work, and the method gap this leaves open

Structural phylogenetics with Foldseek's 3Di alphabet (van Kempen et al., *Nature
Biotechnology* 2022) is an active field, and two papers define its current standard:

- Puente-Lelievre et al. (2024) — 3Di characters as standard phylogenetic characters, with
  IQ-TREE maximum likelihood, a 3Di-specific rate matrix, partitioning, and ultrafast
  bootstrap. Combining amino acids with 3Di best matches a reference structural-distance
  tree and avoids long-branch attraction.
- Fullmer et al. (2025) — 3Di combined with sequence resolves better than either alone, and
  the gain is **weaker in alpha-helical proteins**, because high helical content reduces the
  information 3Di alignments carry.

**Both ask whether 3Di improves phylogenetic resolution. Neither asks whether it changes the
downstream evolutionary events you read off the tree**, which is this repository's question
and the reason it exists.

Two caveats follow, and the second is the important one. Fullmer's alpha-helix result bears
directly on the clades chosen here: MFS transporters are almost entirely helical, so this is
close to the regime where 3Di carries least information. And the method used here — a shared
alignment-free 3-mer cosine distance with neighbour-joining — is deliberately weaker than the
field's, chosen so both alphabets pass through an identical step and no uncalibrated
substitution matrix is compared against another. That control is worth having, but it means
the event-set difference reported above cannot yet be separated from method noise. **Redoing
this with model-based inference and bootstrap support is the next step, and until then the
finding is a motivation rather than a result.**

### More

- [Analysis](ANALYSIS.md) — what was done and why, including two failure modes that would have gone wrong quietly
- [Results](results/RESULTS.md) — full results
- [Design](docs/DESIGN.md) — where the field disagrees with itself, the layout, and the traps this avoids
