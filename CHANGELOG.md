# Changelog
## 2025-04-16
### Changed
- Update Kofam version (`2023-04-01` -> `2025-01-01`), remove threshold scale (`0.75`) for Kofam parsing.
- Omit `env_nr`, use only `nr` for protein database construction.
- Use more strict `diamond` evalue cutoff (`1e-15` -> `1e-25`).
- Adjust marker gene sets (bacteria: `s9` -> `l13`; archaea: `l2` -> `l4e`, `l18e` -> `l44e`, `s28e` -> `s24e`) based on RefSeq re-evaluation.

### Fixed
- Fix NCBI taxonomy (`superkingdom` -> `domain`).