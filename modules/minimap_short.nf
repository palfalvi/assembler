process minimap2_sr {

  label 'assembly'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/bwa", mode: 'copy'

  input:
    tuple val(sample_id), file(reads)
    path assembly

  output:
    path "*.bam", emit: bam
    path "*.bam.bai", emit: baidx

  script:
    """
    minimap2 --secondary=no --MD -L -t $task.cpus -ax map-ont $assembly $reads | samtools sort -@ $task.cpus -O BAM - > ${assembly.simpleName}-${sample_id}.bam
    samtools index ${assembly.simpleName}-${sample_id}.bam
    """
}
