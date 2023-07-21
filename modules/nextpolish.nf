process nextpolish {

  label 'assembly'

  conda "$baseDir/conda-envs/nextpolish-env.yaml"

  // publishDir "${params.outdir}/nextpolish", mode: 'copy'

  input:
    path nextpolish
    path assembly
    path reads

  output:
    path "*nextpolish.fasta", emit: assembly

  script:
    """
    ls  $reads > sgs.fofn

    echo -e "task = best\ngenome = ${assembly}\nsgs_fofn = sgs.fofn" > run.cfg

    nextPolish run.cfg
    """
}
