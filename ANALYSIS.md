# Analysis

What was built, why it was built that way, and two failure modes that would have been
silent.

## The question

Almost every large-scale reconstruction of domain-rearrangement history is built on amino
acid sequence alone. Foldseek's 3Di alphabet turns a fold into a 20-letter string that
existing phylogenetic software consumes unchanged. The question is not whether the trees
differ, since they will, but whether the inferred **rearrangement events** differ. A tree
can be reshuffled substantially while every inferred gain, loss and duplication stays the
same, and that outcome would mean 3Di changes the phylogeny without changing the
evolutionary story.

## Design decisions, and the reasoning

**Both alphabets go through an identical method.** The same alignment-free 3-mer cosine
distance, the same neighbour joining. A substitution-matrix alignment would need a 3Di
matrix and an amino acid matrix calibrated against each other, and they are not. Any
difference in the resulting trees would then be partly a difference between two matrices.
Sharing one method removes that confound, at a cost in resolution which is stated rather
than hidden.

Both alphabets have twenty letters, so the feature spaces are the same size and neither is
advantaged by the representation.

**A globular control clade.** MFS transporters are where rearrangement is claimed to be most
active. Trypsin-domain proteases are the control, where a change would be harder to
attribute to biology.

**Events are read off cherries.** A sister pair is the one place two extant architectures
can be compared without reconstructing an ancestor. Restricting to cherries keeps the
comparison free of an ancestral-state model whose assumptions would otherwise be doing part
of the work.

**Events have no direction.** Two sister paralogues have no ancestor-descendant order
without an outgroup, so a gain in one reading is a loss in the other. Changes are classed
as terminal indel, internal indel or complex, and never as an addition or a deletion.

**Architectures are ordered, not sets.** A terminal indel and an internal indel are
different events and a set cannot tell them apart, so domains are sorted by start position.

**pLDDT masking is a first-class step.** Where AlphaFold is unconfident the backbone is a
guess, so the 3Di letter is a guess about a guess. Disordered linkers, the regions this
question most wants to discuss, are exactly the low-confidence ones. Masking is not a
post-hoc robustness check.

## Two failure modes that would have been silent

**AlphaFold model versions move.** The `v4` URL template that was current when this repo was
designed now returns `NoSuchKey`. A hardcoded template fills the cache with 127-byte XML
error documents, and Foldseek then fails on them for a reason that looks unrelated to the
actual problem. URLs are resolved through the API and a payload under 1 kB is rejected as
not-a-model.

**A masked residue is not a letter.** Any *k*-mer window touching an `X` is dropped rather
than counted. Counting `X`-containing k-mers turns disordered regions into their own signal,
which is precisely the artefact masking exists to remove. A test pins it.

A third turned up while testing the CLI. **A missing Foldseek used to be discovered once per
protein**, every protein was skipped for the same reason, and the run died complaining about
having too few proteins for a tree. True, and silent about the cause. There is now a
preflight check.

## What was measured

58 human proteins (29 transporter, 29 globular), 20 distinct architectures.

- **Robinson-Foulds 80, normalised 0.727.** Most of the non-trivial splits in one tree are
  absent from the other. There is no bootstrap or resampling baseline for this number.
- **Sister pairs.** Only 9 are shared between the two trees (Jaccard 0.30).
- **Architecture changes across sister pairs.** 11 of 21 in the sequence tree and 5 of 18
  in the 3Di tree. The difference is not significant (two-sided Fisher exact p = 0.19).

So the characters change which proteins end up as sisters. This run does not show that
they change how many rearrangements are inferred.

## What is not established

- Whether the event difference is real. It is within what chance allows at this sample
  size, and the complex class includes Pfam family swaps within the MFS clan (CL0015, for
  example PF07690 and PF00083), which are label differences rather than rearrangements.
- How much of the topology difference is noise. Two coarse NJ trees can differ from
  resampled characters alone, and that baseline was not computed.
- Rates. Around twenty cherries per tree cannot estimate how often each event occurs.
- Anything phylogenetic. These are paralogues within one species, so the result is about how
  characters change an inference, not about the evolution of these families.

## What would change the conclusion

Model-based tree inference with the published 3Di substitution matrix, against a standard
amino acid model, with bootstrap support. That would trade the shared-method guarantee for
resolution, and the two runs together would separate "the characters differ" from "the
models differ". Collapsing Pfam families to clans would remove the family-swap events.
