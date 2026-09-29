# domain-arch-3di design notes

Proteins evolve by rearranging modular domains, through terminal additions and deletions,
internal duplications in repeat families, fusion and fission, and, particularly in
transporter families, recombination that produces genuinely novel multidomain
architectures.

**The question.** When the characters change from amino acids to 3Di, do the inferred
*rearrangement events* change, or only the tree topology?

| | |
| --- | --- |
| **Clades** | MFS transporters (IPR020846) plus trypsin-domain serine proteases (IPR001254) as a globular control |
| **Architectures** | Pfam domains from the live InterPro API, ordered by start position |
| **Structures** | AlphaFold DB, with per-residue pLDDT retained |
| **Characters** | Amino acid and 3Di, each through the same 3-mer cosine distance |
| **Readout** | Undirected architecture change across each sister pair, per tree, and the overlap of the two sets of sister pairs |

## Why now, and where the field disagrees with itself

- **Foldtree** (Moi et al., *Nat Struct Mol Biol* 2025,
  [doi:10.1038/s41594-025-01649-8](https://doi.org/10.1038/s41594-025-01649-8)) showed
  structure-derived trees outperform sequence past the twilight zone, and used structural
  phylogenetics to resolve gram-positive communication systems.
- **MBE 2025 (`msaf149`)** concluded that structure-based methods do *not* outperform
  standard sequence methods for large-scale phylogenomics.
- A **general 3Di substitution matrix** (MBE 2025, `msaf124`) and **BEAST 2 support** (Dec
  2025) removed the tooling excuse.
- **"Know Your Alphabet"** (bioRxiv 2026) measured topological variance of 3Di across NMR
  ensembles. The alphabet itself is conformationally noisy.

## Traps the pipeline handles

- **Low-pLDDT regions produce meaningless 3Di.** Disordered linkers are exactly where
  eukaryotic proteins are claimed to expand fastest, and exactly where AlphaFold is least
  confident. Any signal found there is an artefact until it survives pLDDT masking, so
  masking is a first-class step with the threshold recorded in `findings.json`, not a
  post-hoc robustness check.
- **Sister paralogues have no polarity.** Without an outgroup an addition in one reading is
  a deletion in the other, so pair events are undirected.
- **AlphaFold model versions move.** The `v4` URL template that was current when this repo
  was designed now returns `NoSuchKey`, and a hardcoded template fills the cache with
  127-byte XML error documents that Foldseek then fails on for an unrelated-looking reason.
  Download URLs are resolved through the API, and a payload under 1 kB is rejected as
  not-a-model.
- **A masked residue is not a letter.** pLDDT masking replaces low-confidence residues with
  `X`. Any *k*-mer window touching one is dropped rather than counted, because counting
  `X`-containing k-mers turns disordered regions into their own signal, the exact artefact
  the masking exists to remove. A test pins it.

Proteins whose AlphaFold model length disagrees with their UniProt sequence are excluded
rather than aligned by position, and proteins below 50% confident residues are excluded as
a first-class step with the exclusions reported.

## Not done

These were planned and are not implemented.

- **Pinned database versions.** InterPro and UniProt are queried live and their releases
  are not recorded. Pfam and InterPro disagree on boundaries, so an architecture is only
  defined relative to a release. New runs record the accession list and the Foldseek
  version in `findings.json`.
- **Pfam clans.** Families are compared as they are called, so a swap between two families
  of one clan (for example PF07690 and PF00083 in the MFS clan CL0015) counts as a complex
  change.
- **A partitioned amino acid plus 3Di character set.**
- **Gene-model quality filtering.** Truncated or fused gene models create fake terminal
  deletions and fake fusions. The proteins here are reviewed human UniProt entries, and no
  BUSCO or equivalent filter is applied.
- **Template leakage.** Where a structure is a template-based prediction of a close
  homologue, its 3Di adds little information while looking like a second character set.
  This is not checked.
- **Uncertainty.** No bootstrap for either tree and no null model for the event counts.

## Layout

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

No test needs a network or a structure.
