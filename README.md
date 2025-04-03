# ASM Workflow

### About
Nextflow workflow for allele-specific methylation extraction using WBGS Data. Optimized for parallel compute using slurm infrastructure for alignment and leveraging dockerized containers.
Briefly:
* Paired-end FASTQs are adapter trimmed using Trim Galore.
* Paired-end Alignment using Bismark.
* Splitting BAMs by chromosome and adding readgroups for variant calling
* Using BisSNP for variant calling per-chromosome and gold standard filtering from the previous [workflow](https://github.com/TyckoLab/CloudASM)
* Using BAMs and VCFs with [DAMEFinder](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00346-8) for getting allele specific methylation calls. 
* Merging per-chromosome results and annotating with GENE_SYMBOLS and producing a final csv file for downstream analysis.

## Reference Containers
All reference containers referenced in the main file can be found hosted in this dropbox directory [here](https://www.dropbox.com/scl/fo/pcvjjm8jegt0au4amtkek/AGti41sCazjQlzgLohNvlFw?rlkey=zvrt76g43vf56yvc1ymn8s1ks&st=g6aigisq&dl=0). Please change the paths to these containers in the appropriate sections in the main.nf file. 
