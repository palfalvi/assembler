process bwa_index {

  label 'small_job'

  conda "$baseDir/conda-envs/bwa-samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path assembly
    val options

  output:
    path "${assembly}.*", emit: index

  script:
    """
    bwa index $options ${assembly}
    """
}
