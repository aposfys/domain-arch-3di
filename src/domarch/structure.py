"""AlphaFold retrieval, 3Di encoding, and the pLDDT mask the whole result depends on.

3Di describes each residue's local structural environment. Where AlphaFold is unconfident
the backbone geometry is a guess, so the 3Di letter is a guess about a guess. Disordered
linkers -- the regions this repo most wants to say something about -- are precisely the
low-confidence regions. Masking is therefore a first-class step, not a robustness check.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/files"

#: Residues below this pLDDT are treated as unresolved. 70 is the conventional boundary
#: between "confident" and "low"; the value used is written into results/findings.json
#: because every downstream number depends on it.
PLDDT_CONFIDENT = 70.0

#: Foldseek's 20-state structural alphabet.
THREE_DI_ALPHABET = "ACDEFGHIKLMNPQRSTVWY"

#: Character substituted for a masked residue. Chosen to be outside the 3Di alphabet so
#: that a masked position can never be silently scored as a structural match.
MASK_CHARACTER = "X"


def mask_by_plddt(
    three_di: str, plddt: Sequence[float], *, threshold: float = PLDDT_CONFIDENT
) -> str:
    """Replace 3Di letters at low-confidence positions with :data:`MASK_CHARACTER`."""
    if len(three_di) != len(plddt):
        raise ValueError(
            f"length mismatch: {len(three_di)} 3Di characters vs {len(plddt)} pLDDT values"
        )
    return "".join(
        MASK_CHARACTER if score < threshold else char
        for char, score in zip(three_di, plddt, strict=True)
    )


def confident_fraction(plddt: Sequence[float], *, threshold: float = PLDDT_CONFIDENT) -> float:
    """Fraction of residues at or above the confidence threshold."""
    if not plddt:
        raise ValueError("cannot summarise an empty pLDDT array")
    return sum(1 for score in plddt if score >= threshold) / len(plddt)


def fetch_structure(uniprot_id: str, out_dir: Path) -> Path:
    """Download one AlphaFold model, cached. pLDDT is read from the B-factor column."""
    raise NotImplementedError("milestone 1: AlphaFold retrieval")


def encode_3di(structure_path: Path) -> tuple[str, list[float]]:
    """Return the 3Di string and per-residue pLDDT for one structure."""
    raise NotImplementedError("milestone 1: Foldseek 3Di encoding")
