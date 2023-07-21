process minimap2_genomes {

  label 'assembly'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/gala/preliminary_comparison", mode: 'copy'

  input:
    path genome

  output:
    path "*.paf", emit: map

  script:
    """
    minimap2 -x asm5 -DP -t ${task.cpus} $genome $genome > ${genome.simpleName}vs${genome.simpleName}.paf
    """
}
