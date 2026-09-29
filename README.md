# domain-arch-3di
Does Foldseek 3Di change the domain rearrangement events read off a tree, or only the tree topology? A pilot on 58 human proteins.

[![CI](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

```
make install
export FOLDSEEK_BIN=/path/to/foldseek
python3 -m domarch.cli fetch --per-clade 30       # UniProt + InterPro
python3 -m domarch.cli analysis --per-clade 30    # AlphaFold, 3Di, both trees, RESULTS.md
make test                                         # no network, no structure
```

### What it does

58 human proteins in two clades, 29 MFS transporters (IPR020846) and 29 trypsin-domain
serine proteases (IPR001254) as a globular control. AlphaFold models, Foldseek 3Di with
pLDDT masking, and the same alignment-free 3-mer cosine distance and neighbour-joining
step for both alphabets. Architecture changes are read off sister pairs (cherries), where
two extant architectures can be compared without reconstructing an ancestor.

### Result so far

The two trees differ. Robinson-Foulds is 80 (normalised 0.727), and only 9 sister pairs
are shared (Jaccard 0.30). There is no bootstrap or resampling baseline, so how much of
that difference two coarse NJ trees would show from noise alone is not known.

| | Sequence tree | 3Di tree |
| --- | ---: | ---: |
| Cherries | 21 | 18 |
| with an architecture change | 11 | 5 |

The gap in events is not significant (two-sided Fisher exact p = 0.19). This run shows
that the characters change which proteins end up as sisters. It does not show that 3Di
implies fewer rearrangements. Full table in [results/RESULTS.md](results/RESULTS.md).

### Limitations

- Pfam splits the MFS clan (CL0015) into several families, so two single-domain MFS
  proteins with different family calls count as a complex change. The complex counts may
  overstate real rearrangements.
- The distance is deliberately coarse. Both alphabets pass through one identical step, so
  no uncalibrated substitution matrix is compared against another. MFS proteins are
  mostly helical, where 3Di carries least information (Fullmer et al. 2025).
- The run in `results/` did not record its accession list or database releases. New runs
  record the accessions and the Foldseek version in `findings.json`.
- The next step is model-based inference with a 3Di substitution matrix and bootstrap
  support, with Pfam calls collapsed to clans.

### More

- [Analysis](ANALYSIS.md), the reasoning behind each design choice and two silent failure modes
- [Results](results/RESULTS.md), generated from `results/findings.json`
- [Prior work](docs/PRIOR_WORK.md), where this sits relative to current 3Di phylogenetics
- [Design](docs/DESIGN.md), the field context, the layout and what is not done
