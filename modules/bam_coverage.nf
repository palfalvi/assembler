process bam_coverage {

  conda "$baseDir/conda-envs/bwa-samtools-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path bam

  output:
    env avg_depth, emit: coverage

  script:
    """
    avg_depth=`samtools depth ${bam}  |  awk '{sum+=\$3} END { print sum/NR}'`
    """
}
