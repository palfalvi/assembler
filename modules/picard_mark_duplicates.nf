process picard_mark_duplicates {

  label "small_job"

  conda "$baseDir/conda-envs/picard-env.yaml"

  input:
    path bam

  output:
    path "*.deduplicated.bam", emit: bam
    path "*.metrics.txt", emit: metrics

  script:

    """
    picard MarkDuplicates \
      --INPUT $bam \
      --OUTPUT ${bam.simpleName}.deduplicated.bam \
      --METRICS_FILE ${bam.simpleName}.metrics.txt \
      --ASSUME_SORTED true \
      --VALIDATION_STRINGENCY LENIENT \
      --REMOVE_DUPLICATES true
    """
}
