"""Protein sets and their domain architectures.

Two clades, and the pairing is the experiment. A transporter family is where domain
rearrangement is claimed to be most active; a globular family is the control, where a
change in the inferred events would be much harder to attribute to biology.

Domain assignments come from InterPro at a pinned database version. An architecture is
only defined relative to a stated source -- Pfam and InterPro disagree, and "the domain
architecture" is not a database-independent fact.
"""

from __future__ import annotations

import json
import time
import urllib.parse
import urllib.request
from dataclasses import dataclass, field
from pathlib import Path

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
INTERPRO_API = "https://www.ebi.ac.uk/interpro/api"

#: The two clades. `(name, InterPro accession, description)`.
CLADES: tuple[tuple[str, str, str], ...] = (
    ("transporter", "IPR020846", "Major facilitator superfamily domain"),
    ("globular", "IPR001254", "Serine protease, trypsin domain"),
)


@dataclass
class Protein:
    """One protein, its sequence, and its domain architecture."""

    accession: str
    name: str
    clade: str
    sequence: str
    #: Ordered Pfam accessions along the chain. The architecture.
    domains: list[str] = field(default_factory=list)

    @property
    def architecture(self) -> str:
        return "-".join(self.domains)


def _get(url: str, attempts: int = 5) -> dict:
    last: Exception | None = None
    for attempt in range(attempts):
        try:
            request = urllib.request.Request(url, headers={"Accept": "application/json"})
            with urllib.request.urlopen(request, timeout=120) as response:
                return json.loads(response.read().decode("utf-8"))
        except Exception as exc:
            last = exc
            time.sleep(2**attempt)
    raise RuntimeError(f"request failed after {attempts} attempts: {url}") from last


def fetch_clade(
    interpro_id: str, clade: str, *, limit: int = 40, reviewed: bool = True
) -> list[Protein]:
    """Reviewed human proteins carrying an InterPro signature."""
    query = f"(xref:interpro-{interpro_id}) AND (organism_id:9606)"
    if reviewed:
        query += " AND (reviewed:true)"
    params = {
        "query": query,
        "fields": "accession,protein_name,sequence,length",
        "format": "json",
        "size": str(min(limit, 500)),
    }
    payload = _get(f"{UNIPROT_SEARCH}?{urllib.parse.urlencode(params)}")

    proteins: list[Protein] = []
    for entry in payload.get("results", [])[:limit]:
        description = entry.get("proteinDescription", {})
        recommended = description.get("recommendedName", {}).get("fullName", {})
        proteins.append(
            Protein(
                accession=entry["primaryAccession"],
                name=recommended.get("value", entry["primaryAccession"]),
                clade=clade,
                sequence=entry["sequence"]["value"],
            )
        )
    return proteins


def fetch_architecture(accession: str) -> list[str]:
    """Ordered Pfam domain accessions along one protein.

    Ordered by start position, because an architecture is a sequence of domains and not a
    set: a terminal addition and an internal insertion are different events, and a set
    cannot tell them apart.
    """
    url = f"{INTERPRO_API}/entry/pfam/protein/uniprot/{accession}/?page_size=100"
    try:
        payload = _get(url)
    except RuntimeError:
        # A protein with no Pfam match returns 404. That is an empty architecture, not a
        # failure, and it must not be silently confused with a fetch that went wrong --
        # so the caller sees an empty list and the count of these is reported.
        return []

    placed: list[tuple[int, str]] = []
    for result in payload.get("results", []):
        metadata = result.get("metadata", {})
        pfam = metadata.get("accession")
        for match in result.get("proteins", []):
            for location in match.get("entry_protein_locations", []):
                starts = [fragment["start"] for fragment in location.get("fragments", [])]
                if starts and pfam:
                    placed.append((min(starts), pfam))
    placed.sort()
    return [pfam for _, pfam in placed]


def build_dataset(
    out_path: Path, *, per_clade: int = 40, pinned_note: str = "InterPro/Pfam, live"
) -> list[Protein]:
    """Fetch both clades with their architectures, cached."""
    if out_path.exists():
        stored = json.loads(out_path.read_text())
        return [Protein(**row) for row in stored["proteins"]]

    proteins: list[Protein] = []
    for clade, interpro_id, _ in CLADES:
        found = fetch_clade(interpro_id, clade, limit=per_clade)
        for protein in found:
            protein.domains = fetch_architecture(protein.accession)
        proteins.extend(found)
        print(f"{clade}: {len(found)} proteins", flush=True)

    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        json.dumps(
            {
                "source": pinned_note,
                "clades": [
                    {"clade": clade, "interpro": accession, "description": description}
                    for clade, accession, description in CLADES
                ],
                "proteins": [
                    {
                        "accession": p.accession,
                        "name": p.name,
                        "clade": p.clade,
                        "sequence": p.sequence,
                        "domains": p.domains,
                    }
                    for p in proteins
                ],
            },
            indent=1,
        )
    )
    return proteins
