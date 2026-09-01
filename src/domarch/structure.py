"""AlphaFold retrieval, 3Di encoding, and the pLDDT mask the whole result depends on.

3Di describes each residue's local structural environment. Where AlphaFold is unconfident
the backbone geometry is a guess, so the 3Di letter is a guess about a guess. Disordered
linkers -- the regions this repo most wants to say something about -- are precisely the
low-confidence regions. Masking is therefore a first-class step, not a robustness check.
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

#: Resolved through the API rather than a URL template, because model versions move.
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{accession}"

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
    """Download one AlphaFold model, cached. pLDDT is read from the B-factor column.

    The download URL is resolved through the API rather than built from a template. Model
    version numbers change -- the v4 URLs that were current when this was designed now
    return NoSuchKey, and a hardcoded template would have filled the cache with 127-byte
    XML error documents that Foldseek then fails on for a reason unrelated to the actual
    problem. The size check below exists for the same reason.
    """
    import json
    import urllib.request

    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / f"{uniprot_id}.pdb"
    if path.exists() and path.stat().st_size > 1000:
        return path

    with urllib.request.urlopen(
        ALPHAFOLD_API.format(accession=uniprot_id), timeout=120
    ) as response:
        entries = json.loads(response.read().decode("utf-8"))
    if not entries:
        raise FileNotFoundError(f"AlphaFold has no model for {uniprot_id}")

    with urllib.request.urlopen(entries[0]["pdbUrl"], timeout=180) as response:
        payload = response.read()
    if len(payload) < 1000:
        raise OSError(f"{uniprot_id}: AlphaFold returned {len(payload)} bytes, not a model")
    path.write_bytes(payload)
    return path


def encode_3di(structure_path: Path, *, foldseek: str = "foldseek") -> tuple[str, list[float]]:
    """Return the 3Di string and per-residue pLDDT for one structure.

    The two are returned together because a 3Di character is only as trustworthy as the
    confidence of the residue it was computed from, and handing them back separately
    invites using one without the other.
    """
    import subprocess
    import tempfile

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "descriptor.tsv"
        completed = subprocess.run(
            [foldseek, "structureto3didescriptor", str(structure_path), str(out)],
            capture_output=True,
            text=True,
            check=False,
        )
        if completed.returncode != 0 or not out.exists():
            raise RuntimeError(
                f"foldseek failed on {structure_path.name}: {completed.stderr.strip()[:200]}"
            )
        fields = out.read_text().split("\t")
    if len(fields) < 3:
        raise RuntimeError(f"unexpected foldseek output for {structure_path.name}")
    sequence_3di = fields[2].strip()

    plddt = read_plddt(structure_path)
    if len(plddt) != len(sequence_3di):
        raise ValueError(
            f"{structure_path.name}: {len(sequence_3di)} 3Di characters but "
            f"{len(plddt)} pLDDT values; refusing to mask one with the other"
        )
    return sequence_3di, plddt


def read_plddt(structure_path: Path) -> list[float]:
    """Per-residue pLDDT from the B-factor column of an AlphaFold PDB."""
    values: list[float] = []
    seen: set[tuple[str, int]] = set()
    for line in structure_path.read_text().splitlines():
        if not line.startswith("ATOM"):
            continue
        chain = line[21]
        residue = int(line[22:26])
        if (chain, residue) in seen:
            continue
        seen.add((chain, residue))
        values.append(float(line[60:66]))
    return values
