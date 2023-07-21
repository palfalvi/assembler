process pd_cov {

  label 'long_job'

  conda "$baseDir/conda-envs/purge_dups-env.yaml"

  publishDir "${params.outdir}/purge_dups", mode: 'copy'

  input:
  path mapping

  output:
    path 'PB.base.cov', emit: base_cov
    path 'cutoffs', emit: cutoffs

  script:
    """
    pbcstat $mapping

    calcuts PB.stat > cutoffs 2>calcults.log
    """
}
