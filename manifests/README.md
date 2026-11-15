# Release manifests

- release_files.tsv: 537 processed inputs/results, relative paths, byte sizes and SHA-256.
- input_manifest.tsv / server_result_manifest.tsv: inventories shipped inside data archives.
- figure_map.tsv: 84 output mappings, including the original 83 mappings and Table S7.

Run `python scripts/tools/verify_assets.py` after extracting the companion data archives.
The historical workstation staging scripts are in `scripts/local_maintenance/` and are
not necessary to use the release. Downloaded archives must match `release-assets.json`.
