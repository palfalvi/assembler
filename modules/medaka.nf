process medaka {

  label 'assembly'

  conda "$baseDir/conda-envs/medaka-env.yaml"

  publishDir "${params.outdir}/medaka", mode: 'copy'

  input:
    path fastq
    path assembly
    val model

  output:
    path "medaka_polish/consensus.fasta", emit: assembly

  script:

    """
    medaka_consensus \
    -i $fastq \
    -d $assembly \
    -o medaka_polish \
    -t $task.cpus \
    -m $model
    """
}
