process kat {
  label "small_plus"

  conda "$baseDir/conda-envs/kat-env.yaml"

  publishDir "${params.outdir}/kat_plot", mode: 'copy'

  input:
    path genome
    tuple val(sample_id), file(reads)

  output:
    path "*dist_analysis.json", emit: kat_json
    path "*.png", emit: kat_plot

  script:
    """
    kat comp -t ${task.cpus} -o ${genome.simpleName}_kat ${reads} ${genome}
    """
}
