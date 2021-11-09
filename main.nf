
nextflow.enable.dsl = 2

include { FASTP } from './modules/fastp.nf'
include { STAR_INDEX_REFERENCE ; STAR_ALIGN } from './modules/star.nf'
include { SAMTOOLS ; SAMTOOLS_MERGE } from './modules/samtools.nf'
include { CUFFLINKS } from './modules/cufflinks.nf'

log.info """\
         RNAseq differential analysis using NextFlow 
         =============================
         genome: ${params.reference_genome}
         annot : ${params.reference_annotation}
         reads : ${params.reads}
         outdir: ${params.outdir}
         """
         .stripIndent()
 
params.outdir = 'results'


process split_fastq {

    input:
    tuple val(name), path(fastq)

    output:
    tuple val(name), path("${name}_1-${/[0-9]/*params.suffix_length}.fq"), path("${name}_2-${/[0-9]/*params.suffix_length}.fq")

    script:
    """
    cat ${fastq[0]} | split \\
        -a ${params.suffix_length} \\
        -d \\
        -l ${params.num_lines} \\
        - \\
        ${fastq[0].getBaseName()}- \\
        --additional-suffix=".fq"

    cat ${fastq[1]} | split \\
        -a ${params.suffix_length} \\
        -d \\
        -l ${params.num_lines} \\
        - \\
        ${fastq[1].getBaseName()}- \\
        --additional-suffix=".fq" 
    """
}

workflow {

    Channel.fromFilePairs(params.reads) \
	| split_fastq \
	| map { name, fastq, fastq1 -> tuple( groupKey(name, fastq.size()), fastq, fastq1 ) } \
        | transpose() \
        | view()
        | set{ read_pairs_ch }

    FASTP( read_pairs_ch )
    STAR_INDEX_REFERENCE( params.reference_genome, params.reference_annotation )
    STAR_ALIGN( FASTP.out.sample_trimmed, STAR_INDEX_REFERENCE.out, params.reference_annotation )
    SAMTOOLS( STAR_ALIGN.out.sample_sam )
    SAMTOOLS_MERGE( SAMTOOLS.out.sample_bam.collect() )
    CUFFLINKS( SAMTOOLS_MERGE.out.gathered_bam, params.reference_annotation )
    
}
