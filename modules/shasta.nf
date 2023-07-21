process shasta {
  label 'assembly'

  // conda "$baseDir/conda-envs/canu-env.yaml"

  publishDir "${params.outdir}/shasta", mode: 'copy'

  input:
    path fastq

  output:
    path "*.fasta", emit: assembly
    path "*.gfa", emit: gfa

  script:
    """
    #shasta
    """
}
