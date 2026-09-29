.PHONY: install data analysis report test clean clean-data all

PYTHON ?= python3
PER_CLADE ?= 30

all: analysis

## Install the package plus dev tooling. Foldseek is a separate binary, found on PATH or
## through FOLDSEEK_BIN.
install:
	$(PYTHON) -m pip install -e ".[dev]"

## UniProt sequences and InterPro domain architectures for both clades
data:
	$(PYTHON) -m domarch.cli fetch --per-clade $(PER_CLADE)

## AlphaFold models, 3Di, both trees, the event comparison and RESULTS.md
analysis: data
	$(PYTHON) -m domarch.cli analysis --per-clade $(PER_CLADE)

## Re-render RESULTS.md from results/findings.json
report:
	$(PYTHON) -m domarch.cli report

test:
	$(PYTHON) -m pytest -q

clean:
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete the cached dataset and structures
clean-data: clean
	rm -rf data/dataset.json data/structures
