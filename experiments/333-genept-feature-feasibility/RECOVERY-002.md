# 333 recovery 002 — canonical exact-symbol subset, before coverage

September 9, 2026. The original all-uppercase dictionary gate remains
NO-GO. Recovery 001 only corrected the archive filename; it did not change
the schema gate.

`work/name-audit.json` found 133,736 string keys, zero empty keys, 11,294
mixed-case keys, 22,635 whitespace-bearing keys and **10,970 collisions
under uppercasing**. The mapping therefore includes aliases, not only
canonical uppercase symbols. Blind uppercasing or stripping would merge
distinct keys and is explicitly forbidden.

At this amendment, no coverage count or biological model result using these
features has been measured. The source archive and download hashes are
unchanged.

## New, narrowly scoped conversion

Select keys by **exact, case-sensitive string equality** to the public
canonical `data/challenge/gene_names.csv` list plus its target list.
Preserve the canonical spelling, including any mixed-case `orf` symbols.
Do not uppercase, strip, resolve an alias, inspect descriptions to choose a
match, average collisions or use response outcomes. Missing exact symbols
stay missing. There is only one dictionary entry for an exact key.

Validate the **selected subset only** as nonempty strings and fixed-length,
finite, nonzero numeric vectors; convert to float32 and recheck finiteness.
Unused alias entries are not certified by this conversion. This narrows the
validated artifact to the actual response/target vocabulary, rather than
claiming the full release has canonical keys.

Save separate `work/canonical-embeddings.npz` and
`work/canonical-results.json`. Original failure records and the initial
code/manifest remain unchanged. Freeze this amendment, converter, its tests,
name-audit and original download/metadata bindings before conversion.

The original coverage gates are **unchanged**: >=285 of 300 targets and
>=70% of the 18,533 canonical genes. No new model, threshold, response
evaluation, API call or statistical selection. Same 8 GiB memory cap and
serial launch gate. If exact matching fails coverage, do not add alias
rescue based on that result.
