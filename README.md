# domain-arch-3di
Does Foldseek 3Di structural encoding change inferred protein domain rearrangement events, or only tree topology?

[![CI](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml/badge.svg)](https://github.com/aposfys/domain-arch-3di/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)

> **Status: skeleton.** The question below is the hypothesis under test, not a result.

Proteins evolve by rearranging modular domains, and almost every large-scale reconstruction of that history is built on amino acid sequence alone. Two things now make the alternative testable: AlphaFold provides a predicted structure for essentially every known protein, and Foldseek's 3Di alphabet encodes local structural environment as 20 letters — so a fold becomes a string that existing phylogenetic software can consume unchanged.

Nobody has asked the domain-architecture version of the question, which is the version that matters if you want to talk about *how proteins gained their parts* rather than *which protein is whose cousin*.

### Running it
```
make install && make data && make analysis && make test
```

### Layout
```
src/domarch/
  architecture.py   domain architecture parsing and rearrangement event classification
  structure.py      AlphaFold retrieval, 3Di encoding, pLDDT masking
  cli.py            `python -m domarch.cli`
```
Planned: `trees.py` (alignment, partitioned models, tree inference), `compare.py` (event-level comparison between sequence and structure trees).

### Design notes
[Where the field disagrees with itself, and the traps the pipeline is built to avoid](docs/DESIGN.md)
