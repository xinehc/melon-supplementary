# Changelog
## April 2026
### Changed
- Update Kofam version to `2026-01-01`.
- Update Marker gene set (`mg.json`) based on latest RefSeq re-evaluation.
- Keep bash scripts for all.
- Remove RefSeq database.


## April 2025
### Changed
- Update Kofam version (`2023-04-01` -> `2025-01-01`), remove threshold scale (`0.75`) for Kofam parsing.
- Omit `env_nr`, use only `nr` for protein database construction.
- Use more strict `diamond` evalue cutoff (`1e-15` -> `1e-25`).
- Adjust marker gene sets (bacteria: `s9` -> `l13`; archaea: `l2` -> `l4e`, `l18e` -> `l44e`, `s28e` -> `s24e`) based on RefSeq re-evaluation.

### Fixed
- Fix NCBI taxonomy (`superkingdom` -> `domain`).
