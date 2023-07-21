process miniasm {

  label 'assembly'

  conda "$baseDir/conda-envs/miniasm-env.yaml"

  publishDir "${params.outdir}/miniasm", mode: 'copy'

  input:
    path fastq

  output:
    path "miniasm_assembly.fasta", emit: assembly
    path "miniasm_assembly.gfa", emit: gfa

  script:
    """
    minimap2 -x ava-ont -t $task.cpus $fastq $fastq | gzip -1 > reads.paf.gz

    miniasm -f $fastq reads.paf.gz > miniasm_assembly.gfa

    any2fasta miniasm_assembly.gfa > miniasm_assembly.fasta
    """
}
