"""Command line entry point: ``python -m domarch.cli`` or ``domarch``."""

from __future__ import annotations

import argparse
from pathlib import Path

from domarch import __version__
from domarch.structure import PLDDT_CONFIDENT

CHARACTER_SETS = ("aa", "3di", "partitioned")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="domarch",
        description="Sequence vs structure characters for domain architecture evolution",
    )
    parser.add_argument("--version", action="version", version=f"domarch {__version__}")
    parser.add_argument("--data-dir", type=Path, default=Path("data"))
    parser.add_argument("--results-dir", type=Path, default=Path("results"))
    sub = parser.add_subparsers(dest="command", required=True)

    fetch = sub.add_parser("fetch", help="InterPro architectures and AlphaFold structures")
    fetch.add_argument("--family", required=True, help="InterPro or Pfam family accession")
    fetch.add_argument("--control", help="globular control family accession")

    encode = sub.add_parser("encode", help="3Di encoding with pLDDT masking")
    encode.add_argument("--plddt-threshold", type=float, default=PLDDT_CONFIDENT)

    trees = sub.add_parser("trees", help="infer trees under each character set")
    trees.add_argument(
        "--characters", nargs="*", choices=CHARACTER_SETS, default=["aa", "3di"]
    )

    sub.add_parser("compare", help="event-level comparison between the trees")

    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    raise SystemExit(f"'{args.command}' is not implemented yet; see README milestones")


if __name__ == "__main__":
    raise SystemExit(main())
