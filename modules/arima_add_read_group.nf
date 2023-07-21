process arima_add_read_group {

  label "small_job"

  conda "$baseDir/conda-envs/picard-env.yaml"

  input:
    path bam

  output:
    path "*.readgrouped.bam", emit: bam

  script:

    """
    picard AddOrReplaceReadGroups \
      --INPUT $bam \
      --OUTPUT ${bam.simpleName}.readgrouped.bam \
      -ID ${bam.simpleName} \
      -LB ${bam.simpleName} \
      -SM ${bam.simpleName} \
      -PL ILLUMINA \
      -PU none
    """
}
