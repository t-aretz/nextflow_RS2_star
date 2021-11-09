process CUFFLINKS {
    label 'cufflinks'
    publishDir params.outdir
     
    input:
    tuple val(sample_name), path(sorted_bam)
    path(annotation)
    
    output:
    path('transcripts.gtf'), emit: cufflinks_gtf 
    
    script:
    """
    cufflinks -G ${annotation} ${sorted_bam}  
    """
}
