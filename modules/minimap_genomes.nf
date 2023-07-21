process minimap2_genomes {

  label 'assembly'

  conda "$baseDir/conda-envs/minimap-env.yaml"

  // publishDir "${params.outdir}/gala/preliminary_comparison", mode: 'copy'

  input:
    path genome1
    path genome2

  output:
    path "*.paf", emit: map

  script:
    """
    minimap2 -x asm5 -DP -t ${task.cpus} $genome1 $genome2 > ${genome1.simpleName}vs${genome2.simpleName}.paf
    """
}
