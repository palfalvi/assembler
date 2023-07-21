process arima_stats {

  label "small_job"

  conda "$baseDir/conda-envs/samtools-env.yaml"

  input:
    path bam
    path baidx

  output:
    path "*.stats", emit: stats

  script:

    """
    perl $baseDir/scripts/get_stats.pl $bam > ${bam.simpleName}.stats
    """
}
