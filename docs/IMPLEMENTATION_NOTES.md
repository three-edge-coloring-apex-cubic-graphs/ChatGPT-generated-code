# Implementation notes

## Clean-room scope

The implementation follows Appendix B.2–B.4 of *Three-edge coloring apex cubic
graphs* and the free-homomorphism routines cited there from arXiv:2603.24880.
It was developed without consulting the authors' original source code.

## Data parsing

The parser follows `FORMAT.md` and performs structural validation while reading
records.

- File vertex numbers are 1-based; internal identifiers are 0-based.
- An upper degree value of `0` represents infinity.
- On a configuration vertex row, the second value is the number of adjacent
  vertices in the listed rotation. The prescribed degree is that count plus the
  number of incident digons.
- The supplied configuration files encode free completions. Ring vertices are
  removed when constructing an outer extension, and every incidence between a
  ring vertex and an internal vertex becomes a distinct degree-one outer
  endpoint.
- Four supplied configurations begin with a legacy `0 0` empty-record
  sentinel. It is ignored before the actual configuration record.
- Standalone rules may omit a final zero digon count, as allowed by the format.
  Combined-rule and auxiliary-rule records require an explicit digon count so
  consecutive records can be parsed unambiguously.

## Center-free blocking in `allHomImages`

The initial blocking test in Appendix B.4 is unrooted. The implementation does
not choose an arbitrary target vertex as a center. A high-degree vertex of an
earlier configuration may map to any compatible target vertex.

The center-specific blocking overload is retained for combined rules and
cartwheels, where the pseudocode genuinely distinguishes a center. The
uncentered and centered routines therefore have separate interfaces.

## Lemma B.3 search and multiplicity

The direct Appendix B.4 search is available with `--literal-search`. The default
search applies implementation-level pruning and memoization that preserve its
output multiplicity:

1. A dart-pair branch is rejected early only when the full free-homomorphism
   calculation must be empty: the identification would force a self-reverse
   dart/loop, or the degree ranges of the identified heads or tails are
   disjoint.
2. Repeated recursive states and automorphic dart-pair subproblems can be solved
   once. Their complete result lists are replayed for every occurrence.
3. The checked-dart rule of Algorithm B.4.1 is unchanged.
4. Generated islands are not deduplicated by isomorphism before counting or
   reducibility checking.

Consequently, the reported B.3 counters are counts of generated island
occurrences, not counts of isomorphism classes.

## Semi-reducibility checker

The checker implements the semi-D and semi-C definitions using all boundary
colorings, planar full/half Kempe-chain topologies, and deletable edge sets. Its
accelerations include:

- quotienting boundary colorings by global permutations of the three colors;
- cached planar half-Kempe templates;
- bit-mask representations of Kempe components;
- Gray-code traversal of component-switching subcubes;
- early termination after finding an obstruction to a candidate C-reduction;
- caching a verdict for an identical island while still counting every
  occurrence.

Use `--no-reducibility-cache` to recompute the verdict for every occurrence.

## Deterministic ordering

Configurations are sorted by filename. For a requested B.3 range, the driver
still processes the complete preceding prefix to construct `K_smaller`, then
emits output only for the requested range. This makes process-level sharding
consistent with the sequential order `K001.conf`, …, `K915.conf`.
