nextflow.enable.dsl=2

// ─────────────────────────────────────────────
// PARAMS
params.R1_FASTQ         = null
params.R2_FASTQ         = null
params.outdir           = "./results"
params.ref_genome       = "/mnt/mohammadi/home/kganapathy/cloudasm_data/hg19/reference/"
params.threads_bismark  = 6
params.threads_samtools = 16
params.sample_id        = "sample1"
params.TMP_DIR          = "/mnt/mohammadi/home/kganapathy/cache_pac/tmp/"
params.ALL_VARIANTS     = "/mnt/mohammadi/home/kganapathy/cloudasm_data/hg19/variants/All_20180423.vcf.gz"
params.REF_GENOME       = "/mnt/mohammadi/home/kganapathy/cloudasm_data/hg19/reference/human_g1k_v37.fasta"

// ─────────────────────────────────────────────
// MODULE IMPORT
include { chr1_split; chr1_RG; chr1_varcall; chr1_RDS; chr1_annotateRDS } from './modules/chr1.nf'
include { chr2_split; chr2_RG; chr2_varcall; chr2_RDS; chr2_annotateRDS } from './modules/chr2.nf'
include { chr3_split; chr3_RG; chr3_varcall; chr3_RDS; chr3_annotateRDS } from './modules/chr3.nf'
include { chr4_split; chr4_RG; chr4_varcall; chr4_RDS; chr4_annotateRDS } from './modules/chr4.nf'
include { chr5_split; chr5_RG; chr5_varcall; chr5_RDS; chr5_annotateRDS } from './modules/chr5.nf'
include { chr6_split; chr6_RG; chr6_varcall; chr6_RDS; chr6_annotateRDS } from './modules/chr6.nf'
include { chr7_split; chr7_RG; chr7_varcall; chr7_RDS; chr7_annotateRDS } from './modules/chr7.nf'
include { chr8_split; chr8_RG; chr8_varcall; chr8_RDS; chr8_annotateRDS } from './modules/chr8.nf'
include { chr9_split; chr9_RG; chr9_varcall; chr9_RDS; chr9_annotateRDS } from './modules/chr9.nf'
include { chr10_split; chr10_RG; chr10_varcall; chr10_RDS; chr10_annotateRDS } from './modules/chr10.nf'
include { chr11_split; chr11_RG; chr11_varcall; chr11_RDS; chr11_annotateRDS } from './modules/chr11.nf'
include { chr12_split; chr12_RG; chr12_varcall; chr12_RDS; chr12_annotateRDS } from './modules/chr12.nf'
include { chr13_split; chr13_RG; chr13_varcall; chr13_RDS; chr13_annotateRDS } from './modules/chr13.nf'
include { chr14_split; chr14_RG; chr14_varcall; chr14_RDS; chr14_annotateRDS } from './modules/chr14.nf'
include { chr15_split; chr15_RG; chr15_varcall; chr15_RDS; chr15_annotateRDS } from './modules/chr15.nf'
include { chr16_split; chr16_RG; chr16_varcall; chr16_RDS; chr16_annotateRDS } from './modules/chr16.nf'
include { chr17_split; chr17_RG; chr17_varcall; chr17_RDS; chr17_annotateRDS } from './modules/chr17.nf'
include { chr18_split; chr18_RG; chr18_varcall; chr18_RDS; chr18_annotateRDS } from './modules/chr18.nf'
include { chr19_split; chr19_RG; chr19_varcall; chr19_RDS; chr19_annotateRDS } from './modules/chr19.nf'
include { chr20_split; chr20_RG; chr20_varcall; chr20_RDS; chr20_annotateRDS } from './modules/chr20.nf'
include { chr21_split; chr21_RG; chr21_varcall; chr21_RDS; chr21_annotateRDS } from './modules/chr21.nf'
include { chr22_split; chr22_RG; chr22_varcall; chr22_RDS; chr22_annotateRDS } from './modules/chr22.nf'
include { chrX_split; chrX_RG; chrX_varcall; chrX_RDS; chrX_annotateRDS } from './modules/chrX.nf'


// ─────────────────────────────────────────────
// WORKFLOW
workflow {

    if (!params.R1_FASTQ || !params.R2_FASTQ) {
        error "Both --R1_FASTQ and --R2_FASTQ must be provided"
    }

    fastq_pair = Channel.of([file(params.R1_FASTQ), file(params.R2_FASTQ)])
    trimmed       = trimAdapters(fastq_pair)
    aligned_bam   = align(trimmed)
    sorted_bam    = sortAndIndex(aligned_bam)
    
    // chr1
    chr1_split(sorted_bam)       | set { chr1_bam }
    chr1_RG(chr1_bam)            | set { chr1_rg_bam_bai }
    chr1_varcall(chr1_rg_bam_bai) | set { chr1_vcf }
    chr1_RDS(chr1_rg_bam_bai, chr1_vcf) | set { chr1_rds }
    chr1_annotateRDS(chr1_rds)

    // chr2
    chr2_split(sorted_bam)       | set { chr2_bam }
    chr2_RG(chr2_bam)            | set { chr2_rg_bam_bai }
    chr2_varcall(chr2_rg_bam_bai) | set { chr2_vcf }
    chr2_RDS(chr2_rg_bam_bai, chr2_vcf) | set { chr2_rds }
    chr2_annotateRDS(chr2_rds)

    // chr3
    chr3_split(sorted_bam)       | set { chr3_bam }
    chr3_RG(chr3_bam)            | set { chr3_rg_bam_bai }
    chr3_varcall(chr3_rg_bam_bai) | set { chr3_vcf }
    chr3_RDS(chr3_rg_bam_bai, chr3_vcf) | set { chr3_rds }
    chr3_annotateRDS(chr3_rds)

    // chr4
    chr4_split(sorted_bam)       | set { chr4_bam }
    chr4_RG(chr4_bam)            | set { chr4_rg_bam_bai }
    chr4_varcall(chr4_rg_bam_bai) | set { chr4_vcf }
    chr4_RDS(chr4_rg_bam_bai, chr4_vcf) | set { chr4_rds }
    chr4_annotateRDS(chr4_rds)

    // chr5
    chr5_split(sorted_bam)       | set { chr5_bam }
    chr5_RG(chr5_bam)            | set { chr5_rg_bam_bai }
    chr5_varcall(chr5_rg_bam_bai) | set { chr5_vcf }
    chr5_RDS(chr5_rg_bam_bai, chr5_vcf) | set { chr5_rds }
    chr5_annotateRDS(chr5_rds)

    // chr6
    chr6_split(sorted_bam)       | set { chr6_bam }
    chr6_RG(chr6_bam)            | set { chr6_rg_bam_bai }
    chr6_varcall(chr6_rg_bam_bai) | set { chr6_vcf }
    chr6_RDS(chr6_rg_bam_bai, chr6_vcf) | set { chr6_rds }
    chr6_annotateRDS(chr6_rds)

    // chr7
    chr7_split(sorted_bam)       | set { chr7_bam }
    chr7_RG(chr7_bam)            | set { chr7_rg_bam_bai }
    chr7_varcall(chr7_rg_bam_bai) | set { chr7_vcf }
    chr7_RDS(chr7_rg_bam_bai, chr7_vcf) | set { chr7_rds }
    chr7_annotateRDS(chr7_rds)

    // chr8
    chr8_split(sorted_bam)       | set { chr8_bam }
    chr8_RG(chr8_bam)            | set { chr8_rg_bam_bai }
    chr8_varcall(chr8_rg_bam_bai) | set { chr8_vcf }
    chr8_RDS(chr8_rg_bam_bai, chr8_vcf) | set { chr8_rds }
    chr8_annotateRDS(chr8_rds)

    // chr9
    chr9_split(sorted_bam)       | set { chr9_bam }
    chr9_RG(chr9_bam)            | set { chr9_rg_bam_bai }
    chr9_varcall(chr9_rg_bam_bai) | set { chr9_vcf }
    chr9_RDS(chr9_rg_bam_bai, chr9_vcf) | set { chr9_rds }
    chr9_annotateRDS(chr9_rds)

    // chr10
    chr10_split(sorted_bam)       | set { chr10_bam }
    chr10_RG(chr10_bam)           | set { chr10_rg_bam_bai }
    chr10_varcall(chr10_rg_bam_bai) | set { chr10_vcf }
    chr10_RDS(chr10_rg_bam_bai, chr10_vcf) | set { chr10_rds }
    chr10_annotateRDS(chr10_rds)

    // chr11
    chr11_split(sorted_bam)       | set { chr11_bam }
    chr11_RG(chr11_bam)           | set { chr11_rg_bam_bai }
    chr11_varcall(chr11_rg_bam_bai) | set { chr11_vcf }
    chr11_RDS(chr11_rg_bam_bai, chr11_vcf) | set { chr11_rds }
    chr11_annotateRDS(chr11_rds)

    // chr12
    chr12_split(sorted_bam)       | set { chr12_bam }
    chr12_RG(chr12_bam)           | set { chr12_rg_bam_bai }
    chr12_varcall(chr12_rg_bam_bai) | set { chr12_vcf }
    chr12_RDS(chr12_rg_bam_bai, chr12_vcf) | set { chr12_rds }
    chr12_annotateRDS(chr12_rds)

    // chr13
    chr13_split(sorted_bam)       | set { chr13_bam }
    chr13_RG(chr13_bam)           | set { chr13_rg_bam_bai }
    chr13_varcall(chr13_rg_bam_bai) | set { chr13_vcf }
    chr13_RDS(chr13_rg_bam_bai, chr13_vcf) | set { chr13_rds }
    chr13_annotateRDS(chr13_rds)

    // chr14
    chr14_split(sorted_bam)       | set { chr14_bam }
    chr14_RG(chr14_bam)           | set { chr14_rg_bam_bai }
    chr14_varcall(chr14_rg_bam_bai) | set { chr14_vcf }
    chr14_RDS(chr14_rg_bam_bai, chr14_vcf) | set { chr14_rds }
    chr14_annotateRDS(chr14_rds)

    // chr15
    chr15_split(sorted_bam)       | set { chr15_bam }
    chr15_RG(chr15_bam)           | set { chr15_rg_bam_bai }
    chr15_varcall(chr15_rg_bam_bai) | set { chr15_vcf }
    chr15_RDS(chr15_rg_bam_bai, chr15_vcf) | set { chr15_rds }
    chr15_annotateRDS(chr15_rds)

    // chr16
    chr16_split(sorted_bam)       | set { chr16_bam }
    chr16_RG(chr16_bam)           | set { chr16_rg_bam_bai }
    chr16_varcall(chr16_rg_bam_bai) | set { chr16_vcf }
    chr16_RDS(chr16_rg_bam_bai, chr16_vcf) | set { chr16_rds }
    chr16_annotateRDS(chr16_rds)

    // chr17
    chr17_split(sorted_bam)       | set { chr17_bam }
    chr17_RG(chr17_bam)           | set { chr17_rg_bam_bai }
    chr17_varcall(chr17_rg_bam_bai) | set { chr17_vcf }
    chr17_RDS(chr17_rg_bam_bai, chr17_vcf) | set { chr17_rds }
    chr17_annotateRDS(chr17_rds)

    // chr18
    chr18_split(sorted_bam)       | set { chr18_bam }
    chr18_RG(chr18_bam)           | set { chr18_rg_bam_bai }
    chr18_varcall(chr18_rg_bam_bai) | set { chr18_vcf }
    chr18_RDS(chr18_rg_bam_bai, chr18_vcf) | set { chr18_rds }
    chr18_annotateRDS(chr18_rds)

    // chr19
    chr19_split(sorted_bam)       | set { chr19_bam }
    chr19_RG(chr19_bam)           | set { chr19_rg_bam_bai }
    chr19_varcall(chr19_rg_bam_bai) | set { chr19_vcf }
    chr19_RDS(chr19_rg_bam_bai, chr19_vcf) | set { chr19_rds }
    chr19_annotateRDS(chr19_rds)

    // chr20
    chr20_split(sorted_bam)       | set { chr20_bam }
    chr20_RG(chr20_bam)           | set { chr20_rg_bam_bai }
    chr20_varcall(chr20_rg_bam_bai) | set { chr20_vcf }
    chr20_RDS(chr20_rg_bam_bai, chr20_vcf) | set { chr20_rds }
    chr20_annotateRDS(chr20_rds)

    // chr21
    chr21_split(sorted_bam)       | set { chr21_bam }
    chr21_RG(chr21_bam)           | set { chr21_rg_bam_bai }
    chr21_varcall(chr21_rg_bam_bai) | set { chr21_vcf }
    chr21_RDS(chr21_rg_bam_bai, chr21_vcf) | set { chr21_rds }
    chr21_annotateRDS(chr21_rds)

    // chr22
    chr22_split(sorted_bam)       | set { chr22_bam }
    chr22_RG(chr22_bam)           | set { chr22_rg_bam_bai }
    chr22_varcall(chr22_rg_bam_bai) | set { chr22_vcf }
    chr22_RDS(chr22_rg_bam_bai, chr22_vcf) | set { chr22_rds }
    chr22_annotateRDS(chr22_rds)

    // chrX
    chrX_split(sorted_bam)       | set { chrX_bam }
    chrX_RG(chrX_bam)            | set { chrX_rg_bam_bai }
    chrX_varcall(chrX_rg_bam_bai) | set { chrX_vcf }
    chrX_RDS(chrX_rg_bam_bai, chrX_vcf) | set { chrX_rds }
    chrX_annotateRDS(chrX_rds)
}

// ─────────────────────────────────────────────
// PROCESS DEFINITIONS

process trimAdapters {
    tag "${R1.getBaseName()}"

    input:
    tuple path(R1), path(R2)

    output:
    tuple path("*_val_1.fq.gz"), path("*_val_2.fq.gz")

    container '/mnt/mohammadi/home/kganapathy/tyckolab_asm/asm.sif'

    script:
    return """
    ADAPTER_A="AGATCGGAAGAGCACACGTCTGAAC"
    ADAPTER_A2="AGATCGGAAGAGCGTCGTGTAGGGA"

    trim_galore \\
        -a \$ADAPTER_A \\
        -a2 \$ADAPTER_A2 \\
        --quality 30 \\
        --length 40 \\
        --paired \\
        --retain_unpaired \\
        --fastqc \\
        ${R1} \\
        ${R2} \\
        --output_dir .
    """
}



process align {
    tag { r1.getBaseName() }

    input:
    tuple path(r1), path(r2)

    output:
    path("aligned.bam")

    cpus 64
    memory '256 GB'
    time '14d'
    container '/mnt/mohammadi/home/kganapathy/tyckolab_asm/asm.sif'

    script:
    return """
    REF_GENOME="${params.ref_genome}"

    bismark_nozip \\
        -q \\
        --bowtie2 \\
        \$REF_GENOME \\
        -N 1 \\
        -1 ${r1} \\
        -2 ${r2} \\
        --un \\
        --score_min L,0,-0.2 \\
        --bam \\
        --multicore ${params.threads_bismark} \\
        -o .

    mv *_bismark_bt2_pe.bam aligned.bam
    """
}


process sortAndIndex {
    tag "sort"

    input:
    path("merged.bam")

    output:
    path("merged_sorted.bam")
    path("merged_sorted.bam.bai")

    container '/mnt/stsi/stsi4/ASM/samtools/samtools_v1.9-4-deb_cv1.sif'

    script:
    return """
    set -euxo pipefail

    /usr/bin/samtools sort -@ ${params.threads_samtools} -o merged_sorted.bam merged.bam
    /usr/bin/samtools index -@ ${params.threads_samtools} merged_sorted.bam
    """
}

