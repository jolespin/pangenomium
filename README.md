# Pangenomium
Dereplicate genomes and proteins into pangenomes and orthologs

## Dependencies
* [skani](https://github.com/bluenote-1577/skani)
* [mmseqs2](https://github.com/soedinglab/MMseqs2)

## Citations
The methodology used for dereplicating genomes into pangenomes and proteins into orthologs was initially published in [VEBA 2.0](https://academic.oup.com/nar/article/52/14/e63/7697622) which uses [skani](https://www.nature.com/articles/s41592-023-02018-3) for pairwise ANI and [MMseqs2](https://www.nature.com/articles/nbt.3988) for protein clustering.

* Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL. Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing. Nucleic Acids Res. 2024 Aug 12;52(14):e63. doi: 10.1093/nar/gkae528. 

* Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023 Nov;20(11):1661-1665. doi: 10.1038/s41592-023-02018-3.

* Steinegger, M., S??ding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026???1028 (2017). https://doi.org/10.1038/nbt.3988

---

## Usage
### End-to-end for genome and protein-level pangenome clustering
Designed for end-to-end genome and protein-level pangenome clustering.  This pipeline clusters genomes into pangenomes and then clusters proteins within pangenomes.


#### mode: batch
Designed for working with mixed cellular and viral pangenomes. 

```bash
# Default (Batch)
binning_directory="Analysis/veba_output/binning"
genome_manifest_file="Analysis/misc/genomes_table.tsv"
output_directory="Analysis/pangenomium_output"
compile-genomes-table.py -i ${binning_directory} | cut -f1,3,4,5 > ${genome_manifest_file}
# Input: [organism_type, id_genome, genome_filepath, protein_filepath]
pangenomium end-to-end -i ${genome_manifest_file} -m batch -o ${output_directory}
```

#### mode: veba
Designed to work as direct replacement for VEBA's cluster module. 

```bash
# VEBA
binning_directory="Analysis/veba_output/binning"
genome_manifest_file="Analysis/misc/genomes_table.tsv"
output_directory="Analysis/pangenomium_output"
compile-genomes-table.py -i ${binning_directory} > ${genome_manifest_file}
# Input: [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, cds_filepath, gff_filepath]
pangenomium end-to-end -i ${genome_manifest_file} -m veba -o ${output_directory}
```
---

### Clustering genomes into pangenomes

#### mode: batch
Designed for working with mixed cellular and viral pangenomes. 

```bash
# Default (Batch)
binning_directory="Analysis/veba_output/binning"
genome_manifest_file="Analysis/misc/genomes_table.tsv"
output_directory="Analysis/pangenomium_output"
compile-genomes-table.py -i ${binning_directory} | cut -f1,3,4,5 > ${genome_manifest_file}
# Input: [organism_type, id_genome, genome_filepath, protein_filepath]
pangenomium cluster-genomes -i ${genome_manifest_file} -m batch -o ${output_directory}
```

#### mode: veba
Designed to work as direct replacement for VEBA's cluster module. 

```bash
# VEBA
binning_directory="Analysis/veba_output/binning"
genome_manifest_file="Analysis/misc/genomes_table.tsv"
output_directory="Analysis/pangenomium_output"
compile-genomes-table.py -i ${binning_directory} > ${genome_manifest_file}
# Input: [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, cds_filepath, gff_filepath]
pangenomium cluster-genomes -i ${genome_manifest_file} -m veba -o ${output_directory}
```

#### mode: cellular
Defaults set for prokaryotic and/or eukaryotic genomes

```bash
# Cellular
ls path/to/cellular_genomes/*.fa.gz > cellular_genome_filepaths.list
output_directory="Analysis/pangenomium_output"
# Input: path/to/genome.fa[.gz] on each line
pangenomium cluster-genomes -g ${genome_filepaths} -m cellular -o ${output_directory} -x fa.gz
```

Can also pipe filepaths instead:

```bash
ls path/to/cellular_genomes/*.fa.gz | pangenomium cluster-genomes -g ${genome_filepaths} -m cellular -o ${output_directory} -x fa.gz
```

#### mode: viral
Defaults set for viral genomes

```bash
# Viral
ls path/to/viral_genomes/*.fa.gz > viral_genome_filepaths.list
output_directory="Analysis/pangenomium_output"
# Input: path/to/genome.fa[.gz] on each line
pangenomium cluster-genomes -g ${genome_filepaths} -m viral -o ${output_directory} -x fa.gz
```
---

### Clustering proteins from pangenomes into orthologs

#### mode: batch
Designed for working with mixed cellular and viral pangenomes. 

```bash
# Default (Batch)

# Input_1: $genome_manifest_file [organism_type, id_genome, genome_filepath, protein_filepath]
# Input_2: $pangenome_file [id_genome, id_pangenome]
# Output: $genome_manifest_file [id_genome, id_pangenome, protein_filepath]

pangenome_file=${output_directory}/genome_clusters.tsv
pangenome_manifest_file="Analysis/misc/pangenomes_table.tsv"
compile-pangenomes-table.py -i ${genome_manifest_file} -c ${pangenome_file} > ${pangenome_manifest_file}

# Input: $genome_manifest_file [id_genome, id_pangenome, protein_filepath]
pangenomium cluster-proteins-from-pangenomes -i ${pangenome_manifest_file} -m batch -o ${output_directory}
```

#### mode: veba
Designed to work as direct replacement for VEBA's cluster module. 

```bash
# VEBA
...
# Input_1: $genome_manifest_file [organism_type, id_sample, id_genome, genome_filepath, protein_filepath, cds_filepath, gff_filepath]
# Input_2: $pangenome_file [id_genome, id_pangenome]
pangenomium cluster-proteins-from-pangenomes -i ${genome_manifest_file} -c ${pangenome_file} -m veba -o ${output_directory}
```


### Clustering proteins into orthologs

Designed for working with mixed cellular and viral pangenomes. 

```bash
ls path/to/proteins/*.faa.gz | pangenomium cluster-proteins  -o ${output_directory} -x faa.gz
```
