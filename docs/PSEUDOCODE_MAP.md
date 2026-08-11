# Pseudocode-to-source map

## Appendix B.2

| Pseudocode | Implementation |
|---|---|
| B.2.1 `boundaryCompletions` for pseudo-embeddings | `boundary_completions`, `src/b2.cpp` |
| B.2.2 `boundaryCompletions` for triangulations with digons | `boundary_completions`, `src/b2.cpp` |
| B.2.3 `addBoundaryDarts` | `add_boundary_darts`, `src/b2.cpp` |
| B.2.4 `addBoundaryDartsDirectly` | `add_boundary_darts_directly`, `src/b2.cpp` |
| B.2.5 `identifyNeighbors` | `identify_neighbors`, `src/b2.cpp` |
| B.2.6 `linkIncidenceListEnds` | `link_incidence_list_ends`, `src/b2.cpp` |
| B.2.7 `extendFromCutVertices` | `extend_from_cut_vertices`, `src/b2.cpp` |
| B.2.8 `findCutTuples` | `find_cut_tuples`, `src/b2.cpp` |
| B.2.9 `removeRing` | internal part of `extend_from_cut_vertices`, `src/b2.cpp` |
| B.2.10 `fromVRotations` | `from_vertex_rotations`, `src/b2.cpp` |
| B.2.11 `enforceSingleDigonIncidence` | `enforce_single_digon_incidence`, `src/b2.cpp` |
| B.2.12 `twoDigonsIncidentWithSameVertex` | `two_digons_incident_with_same_vertex`, `src/b2.cpp` |
| B.2.13 free homomorphism + single-digon incidence | `free_homomorphism_and_enforce_single_digon_incidence`, `src/b2.cpp` |

## Appendix B.3

| Pseudocode | Implementation |
|---|---|
| B.3.1 `enumDigonIncidentWheels` | `enum_digon_incident_wheels`, `src/b3.cpp` |
| B.3.2 `enumWheels` | `enum_wheels`, `src/b3.cpp` |
| B.3.3 `generateCartwheel` | `generate_cartwheel`, `src/b3.cpp` |
| B.3.4 `neverApply` | `never_apply`, `src/configurations.cpp` |
| B.3.5 lower digon-charge bound | `lower_bound_of_digon_charge`, `src/b3.cpp` |
| B.3.6 `enumDigons` | `enum_digons`, `src/b3.cpp` |
| B.3.7 `fixInRules` | `fix_in_rules`, `src/b3.cpp` |
| B.3.8 `updateDegreeByRule` | `update_degree_by_rule`, `src/b3.cpp` |
| B.3.9 concrete degrees except tail | `concrete_degree_except_tail`, `src/b3.cpp` |
| B.3.10 `prune` | `prune_cartwheel`, `src/b3.cpp` |
| B.3.11 `fixOutRules` | `fix_out_rules`, `src/b3.cpp` |
| B.3.12 `shouldRefine` | `should_refine`, `src/b3.cpp` |
| B.3.13 `refinement` | `refinement`, `src/b3.cpp` |
| B.3.14 `enumPossibleBadWheels` | `enum_possible_bad_wheels`, `src/b3.cpp` |
| B.3.15 verification of one wheel | `verify_no_bad_cartwheels`, `src/b3.cpp` |
| B.3.16 all center degrees | `verify_lemma_b2`, `src/verify.cpp` |

## Appendix B.4

The default Lemma B.3 driver executes these routines with the
correctness-preserving branch rejection and replay memoization described in
`IMPLEMENTATION_NOTES.md`. It does not deduplicate island occurrences or remove
actual dart-pair multiplicity. `--literal-search` disables the B.4 accelerations
and retains the direct nested search for regression comparison.

| Pseudocode | Implementation |
|---|---|
| B.4.1 `allHomImages` | `all_hom_images`, `src/b4.cpp`; its initial blocking test calls the uncentered overload in `src/configurations.cpp` and does not select a target vertex |
| B.4.2 `getWalks` | `get_walks`, `src/b4.cpp` |
| B.4.3 `isPlanar` | `is_planar`, `src/b4.cpp` |
| B.4.4 `hasSeparatingCycle` | `has_separating_cycle`, `src/b4.cpp` |
| B.4.5 `enumCycles` | `enum_cycles`, `src/b4.cpp` |
| B.4.6 `labelDarts` | `label_darts`, `src/b4.cpp` |
| B.4.7 `numSeparatedVertices` | `num_separated_vertices`, `src/b4.cpp` |
| B.4.8 `makeOuterExtension` | `make_outer_extension`, `src/b4.cpp` |
| B.4.9 `findFourDarts` | `find_four_darts`, `src/b4.cpp` |
| B.4.10 `ensureOuterExtension` | `ensure_outer_extension`, `src/b4.cpp` |
| B.4.11 free completion | `free_completion_from_outer_extension`, `src/b4.cpp` |
| B.4.12 island construction | `island_from_free_completion`, `src/b4.cpp` |
| B.4.13 boundary edge indexing | first phase of `island_from_free_completion` |
| B.4.14 pendant edge indexing | second phase of `island_from_free_completion` |
| B.4.15 other edge indexing | third phase of `island_from_free_completion` |
| B.4.16 `constructIsland` | final phase of `island_from_free_completion` |

## Cited dependencies from arXiv:2603.24880

- A.2.1 rooted homomorphism: `homomorphism`, `src/free_hom.cpp`
- A.3.1 dart quotient: `dart_identification`, `src/free_hom.cpp`
- A.4.3–A.4.9 degree issue resolution: `free_homomorphism` and `resolve_degree_issues`, `src/free_hom.cpp`
- A.6 configuration containment and A.7 blocking: `src/configurations.cpp`
- A.8 rule combinations: `combine_rules`, `src/configurations.cpp`
- A.9.1, A.9.3, A.9.4, A.9.12, A.9.13: `always_apply`, charge routines, pruning routines in `src/configurations.cpp` and `src/b3.cpp`
