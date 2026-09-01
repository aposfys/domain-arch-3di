# Analysis

What was built, why it was built that way, and two failure modes that would have been
silent.

## The question

Almost every large-scale reconstruction of domain-rearrangement history is built on amino
acid sequence alone. Foldseek's 3Di alphabet turns a fold into a 20-letter string that
existing phylogenetic software consumes unchanged. The question is not whether the trees
differ — they will — but whether the inferred **rearrangement events** differ. A tree can be
reshuffled substantially while every inferred gain, loss and duplication stays the same;
that outcome would mean 3Di changes the phylogeny without changing the evolutionary story.

## Design decisions, and the reasoning

**Both alphabets go through an identical method.** The same alignment-free 3-mer cosine
distance, the same neighbour joining. A substitution-matrix alignment would need a 3Di
matrix and an amino acid matrix calibrated against each other, and they are not — any
difference in the resulting trees would then be partly a difference between two matrices.
Sharing one method removes that confound, at a cost in resolution which is stated rather
than hidden.

Both alphabets have twenty letters, so the feature spaces are the same size and neither is
advantaged by the representation.

**A globular control clade.** MFS transporters are where rearrangement is claimed to be most
active; trypsin-domain proteases are the control where a change would be harder to attribute
to biology.

**Events are read off cherries.** A sister pair is the one place two extant architectures
can be compared without reconstructing an ancestor. Restricting to cherries keeps the
comparison free of an ancestral-state model whose assumptions would otherwise be doing part
of the work.

**Architectures are ordered, not sets.** A terminal addition and an internal insertion are
different events and a set cannot tell them apart, so domains are sorted by start position.

**pLDDT masking is a first-class step.** Where AlphaFold is unconfident the backbone is a
guess, so the 3Di letter is a guess about a guess — and disordered linkers, the regions this
question most wants to discuss, are exactly the low-confidence ones. Masking is not a
post-hoc robustness check.

## Two failure modes that would have been silent

**AlphaFold model versions move.** The `v4` URL template that was current when this repo was
designed now returns `NoSuchKey`. A hardcoded template fills the cache with 127-byte XML
error documents, and Foldseek then fails on them for a reason that looks unrelated to the
actual problem. URLs are resolved through the API and a payload under 1 kB is rejected as
not-a-model.

**A masked residue is not a letter.** Any *k*-mer window touching an `X` is dropped rather
than counted. Counting `X`-containing k-mers turns disordered regions into their own signal
— precisely the artefact masking exists to remove. A test pins it.

A third was added after the CLI review: **a missing Foldseek used to be discovered once per
protein**, every protein was skipped for the same reason, and the run died complaining about
having too few proteins for a tree. True, and silent about the cause. There is now a
preflight check.

## What was measured

58 human proteins (29 transporter, 29 globular), 20 distinct architectures.

- **Robinson-Foulds 0.727 normalised.** About three-quarters of the bipartitions in one tree
  are absent from the other.
- **Events: 11 rearrangements from the sequence tree, 5 from the structure tree.** Only 9
  sister pairs are shared (Jaccard 0.30).

So **both change**. It is not the case that 3Di reshuffles the phylogeny while leaving the
evolutionary story intact, and here the structural characters imply less than half as many
rearrangements.

**Consistency check:** of the 9 cherries both trees found, 0 were classified differently.
That has to be zero — the same pair of architectures must classify the same way whichever
tree produced it — and it is reported rather than assumed, because a non-zero value would
mean the classifier depends on something beyond the architectures.

## What is not established

- Rates. Around twenty cherries per tree shows the event sets differ; it cannot estimate how
  often each event occurs.
- Anything phylogenetic. These are paralogues within one species, so the result is about how
  characters change an inference, not about the evolution of these families.
- Anything about model-based inference. These trees are alignment-free and correspondingly
  coarse.

## What would change the conclusion

Model-based tree inference with the published 3Di substitution matrix, against a standard
amino acid model. That would trade the shared-method guarantee for resolution, and the two
runs together would separate "the characters differ" from "the models differ".
