# domain-arch-3di — design notes

Proteins evolve by rearranging modular domains: terminal additions and deletions, internal
duplications in repeat families, fusion and fission, and — particularly in transporter
families — recombination that produces genuinely novel multidomain architectures.

**The question:** when the characters change from amino acids to 3Di, do the inferred
*rearrangement events* change, or only the tree topology?

| | |
| --- | --- |
| **Clade** | One transporter family (MFS or ABC) plus a globular control family |
| **Architectures** | InterPro / Pfam domain assignments |
| **Structures** | AlphaFold DB, with per-residue pLDDT retained |
| **Characters** | Amino acid · 3Di · partitioned (both) |
| **Readout** | Fraction of gain/loss/duplication events that change identity or placement |

## Why now — the field disagrees with itself

- **Foldtree** (*Nature* 2023) showed structure-derived trees outperform sequence past the
  twilight zone, and a 2025 *Nat Struct Mol Biol* study used structural phylogenetics to
  resolve gram-positive communication systems.
- **MBE 2025 (`msaf149`)** concluded that structure-based methods do *not* outperform
  standard sequence methods for large-scale phylogenomics.
- A **general 3Di substitution matrix** (MBE 2025, `msaf124`) and **BEAST 2 support** (Dec
  2025) removed the tooling excuse.
- **"Know Your Alphabet"** (bioRxiv 2026) measured topological variance of 3Di across NMR
  ensembles — the alphabet itself is conformationally noisy.

## Traps this pipeline is built to avoid

- **Low-pLDDT regions produce meaningless 3Di.** Disordered linkers are exactly where
  eukaryotic proteins are claimed to expand fastest, and exactly where AlphaFold is least
  confident. Any signal found there is an artefact until it survives pLDDT masking, so
  masking is a first-class step with the threshold recorded in `findings.json`, not a
  post-hoc robustness check.
- **Domain boundaries are model boundaries, not physical ones.** Pfam and InterPro
  disagree; an architecture is only defined relative to a stated database version, which is
  pinned in the data step.
- **Gene prediction errors masquerade as architecture change.** Truncated or fused gene
  models create fake terminal deletions and fake fusions. Architectures derived from
  proteomes with poor BUSCO completeness are excluded, and the exclusion is reported.
- **A 3Di string is not independent evidence from its own structure.** Where a structure is
  a template-based prediction of a close homologue, its 3Di adds no information while
  looking like a second character set. Predicted structures are checked for template
  leakage before being treated as independent.
