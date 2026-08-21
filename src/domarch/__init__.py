"""Sequence versus structure characters for reconstructing domain architecture evolution.

The repository is arranged so that the *event* layer (`architecture`) is independent of the
*character* layer (`structure`, `trees`). That separation is the experiment: the same
rearrangement classifier consumes trees built from amino acids and trees built from 3Di, so
any difference in the conclusion is attributable to the characters and nothing else.
"""

__version__ = "0.1.0"

__all__ = ["__version__"]
