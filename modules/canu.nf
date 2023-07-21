process canu {

  label 'assembly'

  conda "$baseDir/conda-envs/canu-env.yaml"

  publishDir "${params.outdir}/canu", mode: 'copy'

  input:
    path fastq
    val genome_size

  output:
    path "*.contigs.fasta", emit: assembly
    path "*.unitigs.gfa", emit: gfa

  script:
    """
    canu \
    -p  canu_assembly \
    -d canu_out \
    genomeSize=$genome_size \
    useGrid='false' \
    -nanopore $fastq
    """
}
