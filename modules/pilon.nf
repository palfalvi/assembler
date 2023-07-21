process pilon {
  label "assembly"

  conda "$baseDir/conda-envs/pilon-env.yaml"

  publishDir "${params.outdir}/short_polished", mode: 'copy'

  input:
    path genome
    path bam
    path baidx

  output:
    path "pilon*.fasta", emit: assembly

  script:

    """
    pilon \
    -Xmx${task.memory.toGiga()}g \
    --genome $genome \
    --bam $bam \
    --diploid \
    --output pilon \
    --threads ${task.cpus}
    """
}
