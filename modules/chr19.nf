nextflow.enable.dsl=2

process chr19_split {
    tag "split_chr19"

    input:
    path("merged_sorted.bam")
    path("merged_sorted.bam.bai")

    output:
    path("chr_19.bam")

    container '/mnt/stsi/stsi4/ASM/samtools/samtools_v1.9-4-deb_cv1.sif'

    script:
    return """
    /usr/bin/samtools view -b merged_sorted.bam 19 > chr_19.bam
    """
}

process chr19_RG {
    tag "rg_chr19"

    input:
    path("chr_19.bam")

    output:
    tuple path("rg_added_chr19.bam"), path("rg_added_chr19.bai")

    container '/mnt/stsi/stsi4/ASM/picard/picard_v1.139_cv3.sif'

    script:
    return """
    java -jar /opt/conda/bin/picard.jar AddOrReplaceReadGroups \\
        I=chr_19.bam \\
        O=rg_added_chr19.bam \\
        RGID=chr19 \\
        RGLB=lib1 \\
        RGPL=ILLUMINA \\
        RGPU=unit1 \\
        RGSM=${params.sample_id} \\
        SORT_ORDER=coordinate \\
        CREATE_INDEX=true \\
        VALIDATION_STRINGENCY=LENIENT
    """
}

process chr19_varcall {
    tag "varcall_chr19"

    input:
    tuple path(bam), path(bai)

    output:
    path("chr19_filtered.vcf")

    container '/mnt/mohammadi/home/kganapathy/tyckolab_asm/asm.sif'

    script:
    return """
    TMP_DIR="${params.TMP_DIR}"
    ALL_VARIANTS="${params.ALL_VARIANTS}"
    REF_GENOME="${params.REF_GENOME}"

    java -Xmx10g -Djava.io.tmpdir=\${TMP_DIR} -jar /genomics-packages/BisSNP-0.82.2/BisSNP-0.82.2.jar \\
        -L 19 \\
        -T BisulfiteGenotyper \\
        -R \${REF_GENOME} \\
        -D \${ALL_VARIANTS} \\
        -I ${bam} \\
        -vfn1 chr19_raw.vcf \\
        -mmq 30 \\
        -mbq 0 \\
        -stand_call_conf 20 \\
        -nt 8 \\
        -out_modes EMIT_HET_SNPS_ONLY

    perl /genomics-packages/BisSNP-0.82.2/sortByRefAndCor.pl \\
        --k 1 \\
        --c 2 \\
        --tmp \${TMP_DIR} \\
        chr19_raw.vcf \\
        /mnt/mohammadi/home/kganapathy/cloudasm_data/hg19/reference/human_g1k_v37.fasta.fai > chr19_sorted.vcf

    java -Xmx10g -Djava.io.tmpdir=\${TMP_DIR} -jar /genomics-packages/BisSNP-0.82.2/BisSNP-0.82.2.jar \\
        -L 19 \\
        -R \${REF_GENOME} \\
        -T VCFpostprocess \\
        -oldVcf chr19_sorted.vcf \\
        -newVcf chr19_filtered.vcf \\
        -snpVcf chr19_sorted.vcf \\
        -o chr19_snp_summary.txt \\
        -maxCov 200 \\
        -minSNPinWind 2
    """
}

process chr19_RDS {
    tag "rds_chr19"

    input:
    tuple path(bam), path(bai)
    path(vcf)

    output:
    path("chr19_derASM.rds")

    container '/mnt/mohammadi/home/kganapathy/tyckolab_asm/bioconductor-damefinder_1.18.0--r44hdfd78af_0.sif'

    script:
    return """
    Rscript /mnt/stsi/stsi4/ASM_DAMEFINDER_R/generate_RDS.R \\
        "${params.sample_id}" \\
        "${bam}" \\
        "${vcf}" \\
        "19" \\
        "chr19_derASM.rds"
    """
}

process chr19_annotateRDS {
    tag "annotate_chr19"

    input:
    path(rds)

    output:
    path("chr19_derASM_annotated.csv")

    container '/mnt/mohammadi/home/kganapathy/tyckolab_asm/bioconductor-damefinder_1.18.0--r44hdfd78af_0.sif'

    script:
    return """
    Rscript /mnt/stsi/stsi4/ASM_DAMEFINDER_R/annotate_RDS.R \\
        "${rds}" \\
        "chr19_derASM_annotated.csv"
    """
}

