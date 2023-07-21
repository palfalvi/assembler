process multiqc {

  conda "$baseDir/conda-envs/multiqc-env.yaml"

  publishDir "${params.outdir}", mode: 'copy'

  input:
    path('*')
    path config

  output:
    path "*.html"

  script:
    """
    export LC_ALL=en_US.utf8
    multiqc .
    """
}
