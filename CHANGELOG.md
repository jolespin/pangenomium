# CHANGELOG
* [2026.9.15] - Added `nucmer` (MUMmer4) as alternative genome clustering backend via `--genome_clustering_algorithm nucmer` in `cluster-genomes` and `end-to-end`; includes parallel pairwise execution via `--n_concurrent_nucmer_tasks`, dnadiff-based ANI/AF extraction, optional dot plot generation via `--generate_dotplots` (available for both skani and nucmer backends, default format changed to PDF), and automatic archiving of nucmer intermediate files and dot plot auxiliary files (`.gp`, `.fplot`, `.rplot`) by type (`.tar.gz`); added `--keep_temporary` to `cluster-genomes`, `cluster-proteins`, and `cluster-proteins-from-pangenomes` (temporary files removed by default); nucmer/dnadiff per-pair logs are written to the log file instead of the console; renamed `--n_threads_skani` to `--n_threads_ani` in `end-to-end`
* [2026.9.14] - Fixed column mismatch in `cluster-genomes` when parsing 5-column batch manifests (`[organism_type, id_genome, genome, proteins, cds]`); the `>=5` column branch incorrectly treated batch format as VEBA format, shifting column indices so genome IDs were set to genome filepaths and protein filepaths (`.faa`) were used in place of genome filepaths (`.fa`), causing `end-to-end` to fail during symlink creation
* [2026.5.19] - Fixed genome ID truncation when filenames contain `.` characters in the ID (e.g., `B02P4-2019_bin.13.fa.gz`) during auto-detection in `cluster-genomes` [issue #14](https://github.com/jolespin/pangenomium/issues/14)
* [2026.4.2] - Added support for genome IDs different than filenames via symlink-based remapping in `cluster-genomes` [issue #10](https://github.com/jolespin/pangenomium/issues/10)
* [2026.3.24] - Changed `ani_edgelist.tsv` to `skani-triangle_results.tsv` which is now kept
* [2026.2.13] - Fixed issue where `organism_type` wasn't be autodetected from manifest when `-m veba`
* [2026.2.13] - Changed default output directories to `pangenomium_output/[module]` with the exception of `end-to-end` which is just `pangenomium_output/`
* [2026.2.13] - Added `--separator` and set default to `_`  [issue #6](https://github.com/jolespin/pangenomium/issues/6)
* [2026.2.13] - Changed `--protein_cluster_prefix` default to `SSPC-`
* [2026.2.13] - Changed `--genome_cluster_prefix` default to `SLC-`
* [2026.2.13] - Fixed slight clustering inconsistency between `VEBA` and `Pangenomium` [issue #4](https://github.com/jolespin/pangenomium/issues/4)
* [2026.2.12] - Added test genomic assets
* [2025.12.29] - Initial release

### Pending
* Export CDS for orthologs if provided
* Add module for collating annotations
* `Anvi'o` support
