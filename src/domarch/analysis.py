"""The analysis: build both trees, compare topology, compare events."""

from __future__ import annotations

import json
import os
import shutil
from dataclasses import asdict, dataclass
from pathlib import Path

from domarch.compare import compare_events
from domarch.data import CLADES, Protein, build_dataset
from domarch.structure import (
    PLDDT_CONFIDENT,
    confident_fraction,
    encode_3di,
    fetch_structure,
    mask_by_plddt,
)
from domarch.trees import cherries, distance_matrix, neighbour_joining, robinson_foulds


@dataclass(frozen=True)
class Encoded:
    """One protein with both its character sets and its confidence."""

    protein: Protein
    amino_acids: str
    three_di: str
    confident_fraction: float


def foldseek_binary() -> str:
    return os.environ.get("FOLDSEEK_BIN") or shutil.which("foldseek") or "foldseek"


def check_foldseek(binary: str) -> None:
    """Fail immediately if Foldseek cannot run.

    Without this the missing binary is discovered once per protein, every protein is
    skipped for the same reason, and the run dies complaining about having too few
    proteins for a tree -- which is true, and says nothing about the actual cause.
    """
    import subprocess

    try:
        completed = subprocess.run(
            [binary, "version"], capture_output=True, text=True, timeout=60, check=False
        )
    except (OSError, subprocess.SubprocessError) as exc:
        raise RuntimeError(
            f"foldseek is not runnable as {binary!r}: {exc}\n"
            "Install it, or point FOLDSEEK_BIN at the binary:\n"
            "  curl -L https://mmseqs.com/foldseek/foldseek-osx-universal.tar.gz | tar xz\n"
            "  export FOLDSEEK_BIN=$PWD/foldseek/bin/foldseek"
        ) from exc
    if completed.returncode != 0:
        raise RuntimeError(f"foldseek at {binary!r} exited {completed.returncode}")


def run(
    data_dir: Path,
    results_dir: Path,
    *,
    per_clade: int = 40,
    plddt_threshold: float = PLDDT_CONFIDENT,
    min_confident_fraction: float = 0.5,
) -> dict:
    """Fetch, encode, build both trees, and write ``findings.json``."""
    results_dir.mkdir(parents=True, exist_ok=True)
    proteins = build_dataset(data_dir / "dataset.json", per_clade=per_clade)

    binary = foldseek_binary()
    check_foldseek(binary)
    usable: list[Encoded] = []
    skipped: dict[str, str] = {}
    for protein in proteins:
        try:
            path = fetch_structure(protein.accession, data_dir / "structures")
            sequence_3di, plddt = encode_3di(path, foldseek=binary)
        except Exception as exc:
            skipped[protein.accession] = f"{type(exc).__name__}: {exc}"
            continue

        if len(plddt) != len(protein.sequence):
            # The AlphaFold model and the UniProt sequence must be the same protein. A
            # length mismatch means they are not, and aligning them by position would be
            # comparing two different molecules.
            skipped[protein.accession] = (
                f"length mismatch: {len(protein.sequence)} residues in UniProt, "
                f"{len(plddt)} in the model"
            )
            continue

        confident = confident_fraction(plddt, threshold=plddt_threshold)
        if confident < min_confident_fraction:
            # Disordered proteins produce mostly meaningless 3Di. Excluding them is a
            # first-class step, not a post-hoc robustness check, and the exclusions are
            # reported.
            skipped[protein.accession] = (
                f"only {confident:.0%} of residues above pLDDT {plddt_threshold:.0f}"
            )
            continue

        usable.append(
            Encoded(
                protein=protein,
                amino_acids=protein.sequence,
                three_di=mask_by_plddt(sequence_3di, plddt, threshold=plddt_threshold),
                confident_fraction=confident,
            )
        )
        print(
            f"  {protein.accession} {protein.clade:<12} "
            f"{confident:.0%} confident, architecture {protein.architecture or '(none)'}",
            flush=True,
        )

    if len(usable) < 4:
        reasons = ", ".join(sorted({r.split(":")[0] for r in skipped.values()})) or "none"
        raise RuntimeError(
            f"only {len(usable)} usable proteins; need at least 4 for a tree.\n"
            f"{len(skipped)} were skipped. Reasons seen: {reasons}.\n"
            "See results/findings.json from a previous run, or lower "
            "--min-confident-fraction if the models are genuinely low confidence."
        )

    names = [row.protein.accession for row in usable]
    architectures = {row.protein.accession: row.protein.architecture for row in usable}

    sequence_tree = neighbour_joining(
        distance_matrix(names, [row.amino_acids for row in usable])
    )
    structure_tree = neighbour_joining(
        distance_matrix(names, [row.three_di for row in usable])
    )

    rf_absolute, rf_normalised = robinson_foulds(sequence_tree, structure_tree)
    sequence_cherries = cherries(sequence_tree)
    structure_cherries = cherries(structure_tree)
    events = compare_events(sequence_cherries, structure_cherries, architectures)

    findings = {
        "configuration": {
            "clades": [
                {"clade": clade, "interpro": accession, "description": description}
                for clade, accession, description in CLADES
            ],
            "per_clade_requested": per_clade,
            "plddt_threshold": plddt_threshold,
            "min_confident_fraction": min_confident_fraction,
            "distance": "alignment-free 3-mer cosine, identical for both alphabets",
            "tree": "neighbour joining (Biopython)",
        },
        "dataset": {
            "usable": len(usable),
            "skipped": len(skipped),
            "skipped_reasons": skipped,
            "per_clade": {
                clade: sum(1 for row in usable if row.protein.clade == clade)
                for clade, _, _ in CLADES
            },
            "architectures_empty": sum(1 for value in architectures.values() if not value),
            "distinct_architectures": len(set(architectures.values())),
        },
        "topology": {
            "robinson_foulds": rf_absolute,
            "robinson_foulds_normalised": round(rf_normalised, 4),
        },
        "events": asdict(events) | {"cherry_jaccard": round(events.cherry_jaccard, 4)},
    }
    (results_dir / "findings.json").write_text(json.dumps(findings, indent=1))
    return findings
