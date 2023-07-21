process bam_merge {

  conda "$baseDir/conda-envs/samtools-env.yaml"
  label 'small_plus'
  //conda "$baseDir/conda-envs/samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path bams
    val name

  output:
    path "*.merged.bam", emit: bam
    path "*.merged.bam.bai", emit: baidx

  script:
    """
    samtools merge ${name}.merged.bam $bams

    samtools index ${name}.merged.bam
    """
}
