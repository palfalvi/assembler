process quast {
  label "small_job"

  conda "$baseDir/conda-envs/quast-env.yaml"

  publishDir "${params.outdir}/quast", mode: 'copy'

  input:
    path genomes // Can be a list of genomes

  output:
    path "quast_out/report.tsv", emit: summary
    path "*", emit: quast

  script:
    def ref       = params.quast_reference  ? "-r ${params.quast_reference}"        : ""
    def features  = params.quast_features   ? "--features ${params.quast_features}" : ""

    """
    quast.py \
    $ref \
    $features \
    --large \
    --eukaryote \
    --threads ${task.cpus} \
    -o quast_out \
    --no-plots \
    --no-icarus \
    --no-html \
    $genomes
    """
}
