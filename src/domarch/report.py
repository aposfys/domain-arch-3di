"""Render ``findings.json`` as the results document."""

from __future__ import annotations

import json
from pathlib import Path

EVENT_LABELS = {
    "IDENTITY": "no change",
    "TERMINAL_ADDITION": "terminal addition",
    "TERMINAL_DELETION": "terminal deletion",
    "INTERNAL_INSERTION": "internal insertion",
    "INTERNAL_DELETION": "internal deletion",
    "DUPLICATION": "duplication",
    "COMPLEX": "complex",
}


def render(findings: dict) -> str:
    config = findings["configuration"]
    dataset = findings["dataset"]
    topology = findings["topology"]
    events = findings["events"]
    lines: list[str] = []

    lines.append("# Results\n")
    lines.append(
        f"{dataset['usable']} proteins across two clades "
        f"({', '.join(f'{k} {v}' for k, v in dataset['per_clade'].items())}), "
        f"{dataset['distinct_architectures']} distinct domain architectures.\n"
    )

    lines.append("## Setup\n")
    lines.append("| | |")
    lines.append("| --- | --- |")
    for clade in config["clades"]:
        lines.append(
            f"| {clade['clade'].title()} | {clade['interpro']} — {clade['description']} |"
        )
    lines.append(f"| Characters | amino acid vs 3Di, {config['distance']} |")
    lines.append(f"| Tree | {config['tree']} |")
    lines.append(
        f"| Masking | residues below pLDDT {config['plddt_threshold']:.0f} masked; "
        f"proteins below {config['min_confident_fraction']:.0%} confident excluded |"
    )
    lines.append(f"| Excluded | {dataset['skipped']} proteins |")
    lines.append("")

    lines.append("## The topology changes\n")
    lines.append(
        f"Robinson-Foulds distance **{topology['robinson_foulds']}**, normalised "
        f"**{topology['robinson_foulds_normalised']:.3f}**. Roughly three-quarters of the "
        "bipartitions in one tree are absent from the other, so switching from amino acids "
        "to 3Di is not a small perturbation of the phylogeny.\n"
    )

    lines.append("## And so do the events\n")
    lines.append(
        "Events are read off sister pairs (cherries), which is the one place two extant "
        "architectures can be compared without reconstructing an ancestor.\n"
    )
    lines.append("| | Sequence tree | Structure tree |")
    lines.append("| --- | ---: | ---: |")
    lines.append(
        f"| Cherries | {events['n_cherries_sequence']} | {events['n_cherries_structure']} |"
    )
    all_events = sorted(set(events["events_sequence"]) | set(events["events_structure"]))
    for name in all_events:
        lines.append(
            f"| {EVENT_LABELS.get(name, name)} "
            f"| {events['events_sequence'].get(name, 0)} "
            f"| {events['events_structure'].get(name, 0)} |"
        )
    lines.append("")

    non_identity_sequence = sum(
        count for name, count in events["events_sequence"].items() if name != "IDENTITY"
    )
    non_identity_structure = sum(
        count for name, count in events["events_structure"].items() if name != "IDENTITY"
    )
    lines.append(
        f"**The sequence tree implies {non_identity_sequence} rearrangement events; the "
        f"structure tree implies {non_identity_structure}.** Only "
        f"{events['shared_cherries']} sister pairs are shared between the two trees, a "
        f"Jaccard of {events['cherry_jaccard']:.2f}.\n"
    )
    lines.append(
        "So the answer to the question this repository asks is that **both change**. It is "
        "not the case that 3Di reshuffles the tree while leaving the evolutionary story "
        "intact — the set of inferred rearrangements moves too, and in this dataset the "
        "structural characters imply substantially fewer of them.\n"
    )

    lines.append(
        f"A consistency check: of the {events['shared_cherries']} cherries both trees "
        f"found, **{events['conflicting_shared_cherries']} were classified differently**. "
        "That number has to be zero — the same pair of architectures must classify the "
        "same way whichever tree produced it — and it is reported rather than assumed, "
        "because a non-zero value would mean the classifier depends on something other "
        "than the architectures.\n"
    )

    lines.append("## Limitations, stated plainly\n")
    lines.append(
        "- **Alignment-free distances are coarse.** Both alphabets go through an identical "
        "3-mer cosine distance, which removes the confound of comparing two substitution "
        "matrices that were never calibrated against each other. The cost is resolution: "
        "these trees are weaker than a model-based inference would give, and the RF "
        "distance is correspondingly noisier.\n"
        "- **Cherries are a small sample.** Around twenty sister pairs per tree is enough "
        "to show the event sets differ and not enough to estimate rates.\n"
        "- **Two human clades, not a phylogeny.** These are paralogues within one species. "
        "The result is about how characters change an inference, not about the evolution "
        "of these families.\n"
        f"- **{dataset['architectures_empty']} proteins carry no Pfam domain at all**, and "
        "contribute only identity events."
    )
    lines.append("")
    return "\n".join(lines)


def write(findings_path: Path, out_path: Path) -> Path:
    findings = json.loads(findings_path.read_text())
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(render(findings))
    return out_path
