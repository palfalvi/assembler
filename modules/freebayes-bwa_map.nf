process freebayes_bwa {

  label 'assembly'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx
    val avg_depth, emit: avg_depth

  script:
    """
    # bwa index
    bwa index ${assembly}

    # bwa map
    bwa mem -t ${task.cpus} ${assembly} ${reads} | samtools sort -@ ${task.cpus} -O BAM - > ${assembly.simpleName}.bam

    samtools index ${assembly.simpleName}.bam

    avg_depth=`samtools depth ${assembly.simpleName}.bam  |  awk '{sum+=\$3} END { print sum/NR}'`
    """
}
