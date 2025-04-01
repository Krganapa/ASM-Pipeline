# ASM Workflow

### About
Nextflow workflow for allele-specific methylation extraction using WBGS Data. Optimized for parallel compute using slurm infrastructure for alignment and leveraging dockerized containers.
Briefly:
* Paired-end FASTQs are adapter trimmed using Trim Galore.
* Sharded into a number of pieces for quicker alignment
* Aligning each shard using Bismark. 
* Merging all BAMS together 
* Splitting BAMs by chromosome and adding readgroups for variant calling
* Using BisSNP for variant calling per-chromosome and gold standard filtering from the previous [workflow](https://github.com/TyckoLab/CloudASM)
* Using BAMs and VCFs with [DAMEFinder](https://epigeneticsandchromatin.biomedcentral.com/articles/10.1186/s13072-020-00346-8) for getting allele specific methylation calls. 
* Merging per-chromosome results and annotating with GENE_SYMBOLS and producing a final csv file for downstream analysis.

