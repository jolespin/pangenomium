# Pangenomium

Scalable pangenomics toolkit for clustering genomes and proteins across large datasets

> [!CAUTION]
> This project is in developmental stages and designed to reproduce the outputs of [VEBA's clustering module](https://github.com/jolespin/veba/blob/main/bin/cluster.py).
> Documentation on GitHub may not be up-to-date during development.

## Installation

```bash
mamba create -n pangenomium -c conda-forge -c bioconda skani mmseqs2 mummer4 gnuplot setuptools 'python>=3.9' -y
mamba activate pangenomium
pip install pangenomium
```

## Dependencies
* [skani](https://github.com/bluenote-1577/skani)
* [mmseqs2](https://github.com/soedinglab/MMseqs2)
* [mummer4](https://github.com/mummer4/mummer)
* [gnuplot](http://www.gnuplot.info/) (required by mummerplot for dot plot rendering)

## Citations
The methodology used for dereplicating genomes into pangenomes and proteins into orthologs was initially published in [VEBA 2.0](https://academic.oup.com/nar/article/52/14/e63/7697622) which uses [skani](https://www.nature.com/articles/s41592-023-02018-3) for pairwise ANI and [MMseqs2](https://www.nature.com/articles/nbt.3988) for protein clustering.

* Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL. Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing. Nucleic Acids Res. 2024 Aug 12;52(14):e63. doi: 10.1093/nar/gkae528. 

* Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023 Nov;20(11):1661-1665. doi: 10.1038/s41592-023-02018-3.

* Steinegger, M., Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3988

* Marçais G, Delcher AL, Phillippy AM, Coston R, Salzberg SL, Zimin A. MUMmer4: A fast and versatile genome alignment system. PLoS computational biology. 2018 Jan 26;14(1):e1005944.

## Quick Start
> [!NOTE]
> Recommended approach is to use `--genome_clustering_algorithm skani` for larger datasets and `nucmer` for smaller datasets where precision matters

```bash
# End-to-end workflow
pangenomium end-to-end \
  -i genomes_table.tsv \
  --mode veba \
  -o output_directory \
  --n_threads_ani 8 \
  --n_concurrent_mmseqs_tasks 4

# Or run individual steps
pangenomium cluster-genomes -i genome_list.txt -o genome_clustering
pangenomium cluster-proteins-from-pangenomes \
  -g genomes_table.tsv \
  -c genome_clustering/output/genomes_to_pangenomes.tsv.gz \
  -o protein_clustering
```

## Commands

### cluster-genomes

Cluster genomes by ANI using skani (default) or nucmer.

```bash
# Using skani (default, fast)
pangenomium cluster-genomes \
  -i INPUT \
  -o OUTPUT_DIR \
  --n_threads 8 \
  --ani_threshold 95.0

# Using nucmer (pairwise, with dot plots)
pangenomium cluster-genomes \
  -i INPUT \
  -o OUTPUT_DIR \
  --genome_clustering_algorithm nucmer \
  --n_threads 2 \
  --n_concurrent_nucmer_tasks 4 \
  --ani_threshold 95.0 \
  --generate_dotplots
```

**Input formats:**
- Simple list: One genome path per line
- Batch format: TSV with columns `[organism_type, id_genome, genome_filepath, protein_filepath]`
- VEBA format: TSV with columns `[organism_type, id_sample, id_genome, genome_filepath, protein_filepath, cds_filepath, gff_filepath]`

**Outputs:**
- `genomes_to_pangenomes.tsv.gz` - Genome to pangenome cluster assignments
- `serialization/genome_clusters.graph.pkl.gz` - NetworkX graph of genome relationships
- `serialization/genome_clusters.dict.pkl.gz` - Dictionary mapping genomes to clusters
- `representatives/genome_representatives.tsv.gz` - Representative genome for each cluster
- `dotplots/` - Dot plot images (only when `--generate_dotplots` is used; auxiliary files archived as `dotplot_auxiliary_files.tar.gz`)

Temporary files (decompressed genomes, symlinks) in `tmp/` are removed by default after completion. Use `--keep_temporary` to retain them.

**Intermediate files** (in `intermediate/genome_clustering/`):

The ANI edge list (`ani_edgelist_processed.tsv`) is a headerless 5-column TSV used for clustering:

| Column | Description |
|--------|-------------|
| 1 | Reference genome (filepath for skani, genome ID for nucmer) |
| 2 | Query genome (filepath for skani, genome ID for nucmer) |
| 3 | ANI — Average Nucleotide Identity (%) |
| 4 | AF_ref — Alignment fraction of reference genome (%) |
| 5 | AF_query — Alignment fraction of query genome (%) |

For skani, this is derived from `skani-triangle_results.tsv` (first 5 columns, header stripped). For nucmer, values are extracted from dnadiff reports: ANI from `AvgIdentity` (1-to-1 or M-to-M) and alignment fractions from `AlignedBases` percentages.

When using the nucmer backend, per-pair intermediate files are automatically archived by type in `intermediate/genome_clustering/archives/`:

| Archive | Contents |
|---------|----------|
| `nucmer_results.delta.tar.gz` | `.delta` — raw nucmer alignments |
| `nucmer_results.filtered_delta.tar.gz` | `.1delta`, `.mdelta` — filtered alignments |
| `nucmer_results.reports.tar.gz` | `.report` — dnadiff summary stats |
| `nucmer_results.coords.tar.gz` | `.1coords`, `.mcoords` — alignment coordinates |
| `nucmer_results.snps.tar.gz` | `.snps` — SNP calls |
| `nucmer_results.diff.tar.gz` | `.rdiff`, `.qdiff`, `.unref`, `.unqry` — breakpoints and unaligned regions |

When dot plots are generated, mummerplot auxiliary files are archived in `output/dotplots/`:

| Archive | Contents |
|---------|----------|
| `dotplot_auxiliary_files.tar.gz` | `.gp`, `.fplot`, `.rplot` — gnuplot scripts and plot data files |

### cluster-proteins

Cluster all proteins using MMseqs2.

```bash
pangenomium cluster-proteins \
  -i proteins.faa \
  -o OUTPUT_DIR \
  --n_threads 8 \
  --minimum_identity_threshold 50.0
```

**Outputs:**
- `protein_clusters.tsv.gz` - Protein to orthogroup assignments
- `serialization/protein_clusters.graph.pkl.gz` - NetworkX graph
- `serialization/protein_clusters.dict.pkl.gz` - Dictionary mapping
- `representatives/protein_representatives.tsv.gz` - Representative sequences

### cluster-proteins-from-pangenomes

Cluster proteins within each pangenome independently using MMseqs2.

```bash
pangenomium cluster-proteins-from-pangenomes \
  -g genomes_table.tsv \
  -c genomes_to_pangenomes.tsv.gz \
  -o OUTPUT_DIR \
  --n_threads_per_task 2 \
  --n_concurrent_tasks 4
```

**Alternative input:** Direct pangenome manifest with `-i` instead of `-g` and `-c`.

**Outputs:**
- `output/proteins_to_orthologs.tsv.gz` - All protein to orthogroup assignments
- `output/pangenome_tables/` - Per-pangenome prevalence matrices (genome × orthogroup)

### end-to-end

Complete workflow: genome clustering, protein clustering, and comprehensive output generation.

```bash
pangenomium end-to-end \
  -i genomes_table.tsv \
  --mode veba \
  -o OUTPUT_DIR \
  --n_threads_ani 8 \
  --n_concurrent_mmseqs_tasks 4
```

**Modes:**
- `batch` - Expects 4-column format: `[organism_type, id_genome, genome_filepath, protein_filepath]`
- `veba` - Expects 7-column VEBA format

**Outputs:**

Genome clustering results:
- `genome_clustering/output/genomes_to_pangenomes.tsv.gz`
- `genome_clustering/output/serialization/` - Graph and dictionary objects
- `genome_clustering/output/representatives/` - Representative genomes

Protein clustering results:
- `protein_clustering/output/proteins_to_orthologs.tsv.gz`
- `protein_clustering/output/pangenome_tables/` - Per-pangenome prevalence tables

Comprehensive outputs in `output/`:
- `genome_clusters.tsv.gz` - Genome metadata with cluster assignments
- `protein_clusters.tsv.gz` - Protein metadata with orthogroup assignments
- `identifiers.tsv.gz` - All identifier mappings
- `pangenome_tables/` - Prevalence matrices for each pangenome
- `representatives.faa` - Representative protein sequences
- `core_pangenomes/` - Core pangenome sequences per cluster
- `serialization/` - Graph and dictionary objects

### compile-pangenomes-table

Utility to create pangenome manifest from genome manifest and cluster assignments.

```bash
pangenomium compile-pangenomes-table \
  -i genomes_table.tsv \
  -c genomes_to_pangenomes.tsv.gz \
  -o pangenomes_manifest.tsv
```

## Input Format Examples

### Simple genome list
```
/path/to/genome_001.fasta
/path/to/genome_002.fasta
/path/to/genome_003.fasta
```

### Batch manifest (4 columns)
```
prokaryote	genome_001	/path/to/genome_001.fasta	/path/to/proteins_001.faa
prokaryote	genome_002	/path/to/genome_002.fasta	/path/to/proteins_002.faa
eukaryote	genome_003	/path/to/genome_003.fasta	/path/to/proteins_003.faa
```

### VEBA manifest (7 columns)
```
prokaryote	sample_A	genome_001	/path/to/genome_001.fasta	/path/to/proteins_001.faa	/path/to/cds_001.ffn	/path/to/annotation_001.gff
prokaryote	sample_A	genome_002	/path/to/genome_002.fasta	/path/to/proteins_002.faa	/path/to/cds_002.ffn	/path/to/annotation_002.gff
prokaryote	sample_B	genome_003	/path/to/genome_003.fasta	/path/to/proteins_003.faa	/path/to/cds_003.ffn	/path/to/annotation_003.gff
```

### Pangenome manifest (3 columns)
```
genome_001	PanG-abc123	/path/to/proteins_001.faa
genome_002	PanG-abc123	/path/to/proteins_002.faa
genome_003	PanG-def456	/path/to/proteins_003.faa
```

## Output File Descriptions

### Cluster assignment files
- **genomes_to_pangenomes.tsv.gz** - Two columns: `[id_genome, id_pangenome]`
- **proteins_to_orthologs.tsv.gz** - Two columns: `[id_protein, id_orthogroup]`

### Comprehensive output tables
- **genome_clusters.tsv.gz** - Genome metadata: ID, organism type, cluster, paths, contig count, total length
- **protein_clusters.tsv.gz** - Protein metadata: ID, genome, contig, cluster, length, product
- **identifiers.tsv.gz** - All identifier mappings: genome, sample, pangenome, contig, protein, orthogroup

### Prevalence tables
Per-pangenome matrices in `pangenome_tables/`:
- Rows: Genomes
- Columns: Orthogroups
- Values: Number of proteins from genome in orthogroup

### Representative sequences
- **representatives.faa** - Representative protein for each orthogroup
- **core_pangenomes/*.faa** - Core proteins for each pangenome cluster

### Serialization objects
- **genome_clusters.graph.pkl.gz** - NetworkX graph of genome ANI relationships
- **genome_clusters.dict.pkl.gz** - Dictionary: genome → pangenome cluster
- **protein_clusters.graph.pkl.gz** - NetworkX graph of protein similarity
- **protein_clusters.dict.pkl.gz** - Dictionary: protein → orthogroup

## Algorithm Details

**Genome clustering:**
- Two backends available via `--genome_clustering_algorithm`:
  - `skani` (default): Fast all-vs-all ANI via `skani triangle`. Best for large datasets.
  - `nucmer`: Pairwise ANI via MUMmer4 `nucmer` + `dnadiff`. Better precision for smaller datasets; enables dot plot visualization.
- Both backends produce the same 5-column ANI edge list (see intermediate files above)
- Clusters with `edgelist-to-clusters.py` using connected components (single-linkage)
- Default thresholds: 95% ANI, 50% alignment fraction (relaxed mode: max of ref/query AF must pass)
- Dot plots (`--generate_dotplots`): Available with either backend (default format: PDF). For skani, nucmer is run post-hoc on threshold-passing pairs only. Auxiliary files (`.gp`, `.fplot`, `.rplot`) are archived automatically.

**Protein clustering:**
- Uses MMseqs2 easy-cluster
- Per-pangenome clustering prevents inflation of orthogroup sizes
- Default: 50% identity, 80% coverage (bidirectional)

**Parallelization:**
- Genome clustering (skani): Multi-threaded via `--n_threads`
- Genome clustering (nucmer): Concurrent pairwise comparisons via `--n_concurrent_nucmer_tasks`, each using `--n_threads` threads
- Protein clustering: Multiple concurrent MMseqs2 jobs, each multi-threaded

## Common Workflows

### From VEBA output
```bash
pangenomium end-to-end \
  -i veba_output/genomes_table.tsv \
  --mode veba \
  -o pangenome_analysis \
  --n_threads_ani 16 \
  --n_concurrent_mmseqs_tasks 8 \
  --n_threads_mmseqs_per_task 4
```

### Manual two-step workflow
```bash
# Step 1: Cluster genomes
pangenomium cluster-genomes \
  -i genomes_table.tsv \
  -o step1_genomes \
  --n_threads 16

# Step 2: Cluster proteins per pangenome
pangenomium cluster-proteins-from-pangenomes \
  -g genomes_table.tsv \
  -c step1_genomes/output/genomes_to_pangenomes.tsv.gz \
  -o step2_proteins \
  --n_concurrent_tasks 8 \
  --n_threads_per_task 4
```

### Protein clustering without genome context
```bash
# Concatenate all proteins
cat proteins/*.faa > all_proteins.faa

# Cluster all proteins together
pangenomium cluster-proteins \
  -i all_proteins.faa \
  -o protein_clusters \
  --n_threads 16
```
