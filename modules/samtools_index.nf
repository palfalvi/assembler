process samtools_index {

  label 'small_job'

  conda "$baseDir/conda-envs/samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path bam
    val options

  output:
    path "*.bai", emit: baidx

  script:
    """
    samtools index $options ${bam}
    """
}
