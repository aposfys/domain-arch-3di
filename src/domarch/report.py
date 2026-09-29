"""Render ``findings.json`` as the results document."""

from __future__ import annotations

import json
from pathlib import Path

from domarch.compare import fisher_exact_two_sided

#: Event names as findings files store them, mapped to undirected labels. Sister paralogues
#: have no direction without an outgroup, so the directional names written by earlier runs
#: are merged into the same classes the current classifier emits.
EVENT_LABELS = {
    "IDENTITY": "no change",
    "TERMINAL_INDEL": "terminal indel",
    "TERMINAL_ADDITION": "terminal indel",
    "TERMINAL_DELETION": "terminal indel",
    "INTERNAL_INDEL": "internal indel",
    "INTERNAL_INSERTION": "internal indel",
    "INTERNAL_DELETION": "internal indel",
    "INTERNAL_DUPLICATION": "internal indel",
    "DUPLICATION": "internal indel",
    "COMPLEX": "complex",
}
LABEL_ORDER = ("no change", "terminal indel", "internal indel", "complex")


def merge_events(counts: dict[str, int]) -> dict[str, int]:
    """Event counts keyed by undirected label, in a fixed order."""
    merged: dict[str, int] = dict.fromkeys(LABEL_ORDER, 0)
    for name, count in counts.items():
        label = EVENT_LABELS.get(name, name.lower())
        merged[label] = merged.get(label, 0) + count
    return merged


def render(findings: dict) -> str:
    config = findings["configuration"]
    dataset = findings["dataset"]
    topology = findings["topology"]
    events = findings["events"]
    lines: list[str] = []

    lines.append("# Results\n")
    lines.append(
        f"A pilot run. {dataset['usable']} proteins across two clades "
        f"({', '.join(f'{k} {v}' for k, v in dataset['per_clade'].items())}), "
        f"{dataset['distinct_architectures']} distinct domain architectures.\n"
    )

    lines.append("## Setup\n")
    lines.append("| | |")
    lines.append("| --- | --- |")
    for clade in config["clades"]:
        lines.append(
            f"| {clade['clade'].title()} | {clade['interpro']}, {clade['description']} |"
        )
    lines.append(f"| Characters | amino acid vs 3Di, {config['distance']} |")
    lines.append(f"| Tree | {config['tree']} |")
    lines.append(
        f"| Masking | residues below pLDDT {config['plddt_threshold']:.0f} masked, "
        f"proteins below {config['min_confident_fraction']:.0%} confident excluded |"
    )
    lines.append(f"| Excluded | {dataset['skipped']} proteins |")
    if config.get("foldseek_version"):
        lines.append(f"| Foldseek | {config['foldseek_version']} |")
    lines.append("")

    lines.append("## Topology\n")
    lines.append(
        f"Robinson-Foulds distance **{topology['robinson_foulds']}**, normalised "
        f"**{topology['robinson_foulds_normalised']:.3f}**. Most of the non-trivial splits "
        "in one tree are absent from the other. There is no bootstrap or resampling "
        "baseline, so how much two alignment-free neighbour-joining trees differ from noise "
        "alone is not known.\n"
    )

    lines.append("## Events\n")
    lines.append(
        "Events are read off sister pairs (cherries), the one place two extant "
        "architectures can be compared without reconstructing an ancestor. Sister "
        "paralogues have no direction without an outgroup, so an indel is not split into "
        "an addition and a deletion.\n"
    )
    sequence = merge_events(events["events_sequence"])
    structure = merge_events(events["events_structure"])
    n_sequence = events["n_cherries_sequence"]
    n_structure = events["n_cherries_structure"]
    lines.append("| | Sequence tree | Structure tree |")
    lines.append("| --- | ---: | ---: |")
    lines.append(f"| Cherries | {n_sequence} | {n_structure} |")
    for label in sequence:
        lines.append(f"| {label} | {sequence[label]} | {structure.get(label, 0)} |")
    lines.append("")

    changed_sequence = n_sequence - sequence["no change"]
    changed_structure = n_structure - structure["no change"]
    p_value = fisher_exact_two_sided(
        changed_sequence, sequence["no change"], changed_structure, structure["no change"]
    )
    verdict = (
        "The difference is not significant"
        if p_value >= 0.05
        else "The difference is nominally significant"
    )
    lines.append(
        f"The sequence tree's cherries show {changed_sequence} architecture changes out of "
        f"{n_sequence}, the structure tree's {changed_structure} out of {n_structure}. "
        f"{verdict} (two-sided Fisher exact p = {p_value:.2f}, treating cherries as "
        "independent, which overstates the evidence because the two trees share "
        f"{events['shared_cherries']} of them). Only {events['shared_cherries']} sister "
        f"pairs are shared between the trees, a Jaccard of {events['cherry_jaccard']:.2f}.\n"
    )
    lines.append(
        "So the two alphabets give different trees and largely different sister pairs. "
        "Whether they imply different numbers of rearrangements is not settled by this run.\n"
    )

    lines.append("## Limitations\n")
    lines.append(
        "- **Complex includes Pfam family swaps.** Pfam splits the MFS clan (CL0015) into "
        "several families, for example PF07690 (MFS_1) and PF00083 (Sugar_tr). Two "
        "single-domain MFS proteins with different family calls count as complex although "
        "no domain was gained, lost or moved, so the complex counts may overstate real "
        "rearrangements.\n"
        "- **Alignment-free distances are coarse.** Both alphabets go through an identical "
        "3-mer cosine distance, which removes the confound of comparing two substitution "
        "matrices that were never calibrated against each other. The cost is resolution. "
        "These trees are weaker than a model-based inference would give, and the RF "
        "distance is correspondingly noisier.\n"
        "- **Cherries are a small sample.** Around twenty sister pairs per tree is too few "
        "to estimate rates or to separate a modest difference from chance.\n"
        "- **Two human clades, not a phylogeny.** These are paralogues within one species. "
        "The result is about how characters change an inference, not about the evolution "
        "of these families."
    )
    if not dataset.get("accessions"):
        lines.append(
            "- **Inputs not recorded.** This run did not record its accession list or the "
            "UniProt, InterPro and Foldseek releases. The protein set is the first "
            f"{config['per_clade_requested']} UniProt hits per clade, which depends on the "
            "release, so a rerun can give different numbers."
        )
    if dataset["architectures_empty"]:
        lines.append(
            f"- **{dataset['architectures_empty']} proteins carry no Pfam domain.** A "
            "cherry pairing one of them with a domain-carrying protein counts as a "
            "terminal indel."
        )
    lines.append("")
    return "\n".join(lines)


def write(findings_path: Path, out_path: Path) -> Path:
    findings = json.loads(findings_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render(findings))
    return out_path
