"""Command line entry point: ``python -m domarch.cli`` or ``domarch``."""

from __future__ import annotations

import argparse
from pathlib import Path

from domarch import __version__
from domarch.structure import PLDDT_CONFIDENT


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="domarch",
        description="Sequence vs structure characters for domain architecture evolution",
    )
    parser.add_argument("--version", action="version", version=f"domarch {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser(
        "fetch", help="UniProt sequences and InterPro architectures for both clades"
    )
    fetch.add_argument("--per-clade", type=int, default=30)

    analysis = sub.add_parser(
        "analysis",
        help="structures, 3Di, both trees, and the event-level comparison",
    )
    analysis.add_argument("--per-clade", type=int, default=30)
    analysis.add_argument("--plddt-threshold", type=float, default=PLDDT_CONFIDENT)
    analysis.add_argument(
        "--min-confident-fraction",
        type=float,
        default=0.5,
        help="exclude proteins with less than this fraction of confident residues",
    )

    report = sub.add_parser("report", help="render RESULTS.md from findings.json")

    del report
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    if args.command == "fetch":
        from domarch.data import build_dataset

        proteins = build_dataset(args.data_dir / "dataset.json", per_clade=args.per_clade)
        print(f"{len(proteins)} proteins in {args.data_dir / 'dataset.json'}")
        return 0

    if args.command == "analysis":
        from domarch.analysis import run
        from domarch.report import write

        findings = run(
            args.data_dir,
            args.results_dir,
            per_clade=args.per_clade,
            plddt_threshold=args.plddt_threshold,
            min_confident_fraction=args.min_confident_fraction,
        )
        write(args.results_dir / "findings.json", args.results_dir / "RESULTS.md")
        print(
            f"RF {findings['topology']['robinson_foulds_normalised']:.3f}, "
            f"cherry Jaccard {findings['events']['cherry_jaccard']:.3f}"
        )
        return 0

    if args.command == "report":
        from domarch.report import write

        out = write(args.results_dir / "findings.json", args.results_dir / "RESULTS.md")
        print(f"wrote {out}")
        return 0

    raise SystemExit(f"unknown command {args.command!r}")


if __name__ == "__main__":
    raise SystemExit(main())
