# Manifests / 清单

| File | Content |
|---|---|
| `data_download.tsv` | Archive URL, SHA256 and extraction note for `data/` and `server_result/` |
| `output_manifest.tsv` | Every expected figure and table, with kind and reference size |
| `../data/MANIFEST.tsv` | Input files: role (raw/prepared/precomputed), size, SHA256, source |
| `../server_result/MANIFEST.tsv` | Server results: source, size, SHA256, generating step/parameters |

`output_manifest.tsv` is used by `scripts/tools/check_outputs.py`. To compare a
new run against a reference, pass the reference result directory to
`scripts/tools/compare_results.py`.

`output_manifest.tsv` 供输出检查脚本使用；数值比对使用
`scripts/tools/compare_results.py --ref <reference> --new result`。
