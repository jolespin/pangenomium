# CHANGELOG
* [2026.2.13] - Fixed issue where `organism_type` wasn't be autodetected from manifest when `-m veba`
* [2026.2.13] - Changed default output directories to `pangenomium_output/[module]` with the exception of `end-to-end` which is just `pangenomium_output/`
* [2026.2.13] - Added `--separator` and set default to `_`  [issue #6](https://github.com/jolespin/pangenomium/issues/6)
* [2026.2.13] - Changed `--protein_cluster_prefix` default to `SSPC-`
* [2026.2.13] - Changed `--genome_cluster_prefix` default to `SLC-`
* [2026.2.13] - Fixed slight clustering inconsistency between `VEBA` and `Pangenomium` [issue #4](https://github.com/jolespin/pangenomium/issues/4)
* [2026.2.12] - Added test genomic assets
* [2025.12.29] - Initial release

### Pending
* Add `--protein_cluster_prefix`
* Export CDS for orthologs if provided
* Add `organsim_type` prefix to pangenomes
* Add module for collating annotations
* Add `--separator` argument with default `_` instead of `__`