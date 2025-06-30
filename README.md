# Pangenomium
Dereplicate genomes and proteins into pangenomes and orthologs

## Dependencies
* [skani](https://github.com/bluenote-1577/skani)
* [mmseqs2](https://github.com/soedinglab/MMseqs2)

## Citations
The methodology used for dereplicating genomes into pangenomes and proteins into orthologs was initially published in [VEBA 2.0](https://academic.oup.com/nar/article/52/14/e63/7697622) which uses [skani](https://www.nature.com/articles/s41592-023-02018-3) for pairwise ANI and [MMseqs2](https://www.nature.com/articles/nbt.3988) for protein clustering.

* Espinoza JL, Phillips A, Prentice MB, Tan GS, Kamath PL, Lloyd KG, Dupont CL. Unveiling the microbial realm with VEBA 2.0: a modular bioinformatics suite for end-to-end genome-resolved prokaryotic, (micro)eukaryotic and viral multi-omics from either short- or long-read sequencing. Nucleic Acids Res. 2024 Aug 12;52(14):e63. doi: 10.1093/nar/gkae528. 

* Shaw J, Yu YW. Fast and robust metagenomic sequence comparison through sparse chaining with skani. Nat Methods. 2023 Nov;20(11):1661-1665. doi: 10.1038/s41592-023-02018-3.

* Steinegger, M., Söding, J. MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. Nat Biotechnol 35, 1026–1028 (2017). https://doi.org/10.1038/nbt.3988