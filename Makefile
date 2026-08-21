.PHONY: install data encode trees analysis test clean clean-data all

PYTHON ?= python3
FAMILY ?= PF07690
CONTROL ?= PF00042

all: analysis

## Install the package plus dev tooling. Foldseek comes from environment.yml:
##   conda env create -f environment.yml && conda activate domarch
install:
	$(PYTHON) -m pip install -e ".[dev]"

## InterPro architectures and AlphaFold structures for the family and its control
data:
	$(PYTHON) -m domarch.cli fetch --family $(FAMILY) --control $(CONTROL)

## 3Di encoding, with low-pLDDT residues masked before anything downstream sees them
encode: data
	$(PYTHON) -m domarch.cli encode

## One tree per character set
trees: encode
	$(PYTHON) -m domarch.cli trees --characters aa 3di partitioned

## Event-level comparison: which rearrangements change, and where they sit
analysis: trees
	$(PYTHON) -m domarch.cli compare

test:
	$(PYTHON) -m pytest -q

clean:
	rm -rf results/*
	find . -name __pycache__ -type d -exec rm -rf {} +

## Also delete cached structures and architecture tables
clean-data: clean
	rm -f data/*.cif data/*.tsv data/*.fasta
