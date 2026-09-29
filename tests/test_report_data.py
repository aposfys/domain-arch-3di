"""Results rendering and the InterPro fetch, neither needing a network."""

from __future__ import annotations

import io
import json
import urllib.error
from pathlib import Path

import pytest

from domarch import data, report
from domarch.analysis import check_foldseek

FINDINGS = Path(__file__).resolve().parent.parent / "results" / "findings.json"


def test_committed_findings_render_with_undirected_events():
    text = report.render(json.loads(FINDINGS.read_text()))
    assert "terminal indel | 3 | 1" in text
    assert "internal indel | 1 | 1" in text
    assert "terminal addition" not in text
    assert "internal deletion" not in text
    assert "p = 0.19" in text
    assert "not significant" in text
    assert "—" not in text


def test_empty_architecture_line_only_appears_when_there_are_some():
    findings = json.loads(FINDINGS.read_text())
    assert "carry no Pfam domain" not in report.render(findings)
    findings["dataset"]["architectures_empty"] = 2
    assert "2 proteins carry no Pfam domain" in report.render(findings)


def test_a_protein_with_no_pfam_match_has_an_empty_architecture(monkeypatch):
    def not_found(url, attempts=5):
        raise data.NotFound(url)

    monkeypatch.setattr(data, "_get", not_found)
    assert data.fetch_architecture("P00000") == []


def test_a_failed_fetch_is_not_recorded_as_no_domains(monkeypatch):
    def unreachable(url, attempts=5):
        raise RuntimeError(f"request failed after {attempts} attempts: {url}")

    monkeypatch.setattr(data, "_get", unreachable)
    with pytest.raises(RuntimeError, match="request failed"):
        data.fetch_architecture("P00000")


def test_get_raises_not_found_on_404_without_retrying(monkeypatch):
    calls = []

    def urlopen(request, timeout):
        calls.append(request)
        raise urllib.error.HTTPError(request.full_url, 404, "Not Found", {}, io.BytesIO())

    monkeypatch.setattr(data.urllib.request, "urlopen", urlopen)
    monkeypatch.setattr(data.time, "sleep", lambda seconds: None)
    with pytest.raises(data.NotFound):
        data._get("https://example.invalid/x")
    assert len(calls) == 1


def test_get_reads_an_empty_body_as_an_empty_payload(monkeypatch):
    class Response(io.BytesIO):
        def __enter__(self):
            return self

        def __exit__(self, *args):
            return False

    monkeypatch.setattr(data.urllib.request, "urlopen", lambda request, timeout: Response())
    assert data._get("https://example.invalid/x") == {}


def test_foldseek_preflight_returns_the_version(tmp_path):
    fake = tmp_path / "foldseek"
    fake.write_text("#!/bin/sh\necho 10.941cd33\n")
    fake.chmod(0o755)
    assert check_foldseek(str(fake)) == "10.941cd33"
