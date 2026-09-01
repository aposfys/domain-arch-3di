"""The CLI surface, and the Foldseek preflight.

Without the preflight a missing binary was discovered once per protein, every protein was
skipped for the same reason, and the run died complaining about having too few proteins --
which was true and said nothing about the cause.
"""

from __future__ import annotations

import pytest

from domarch.analysis import check_foldseek
from domarch.cli import build_parser, main


def test_missing_foldseek_is_named_directly():
    with pytest.raises(RuntimeError) as excinfo:
        check_foldseek("/definitely/not/a/binary")
    message = str(excinfo.value)
    assert "foldseek is not runnable" in message
    assert "FOLDSEEK_BIN" in message


def test_report_without_findings_says_what_to_run(tmp_path):
    with pytest.raises(SystemExit) as excinfo:
        main(["--results-dir", str(tmp_path), "report"])
    assert "domarch analysis" in str(excinfo.value)


def test_subcommands_parse():
    parser = build_parser()
    assert parser.parse_args(["fetch"]).command == "fetch"
    assert parser.parse_args(["analysis"]).command == "analysis"
    assert parser.parse_args(["report"]).command == "report"
