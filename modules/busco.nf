process busco {
  tag "$genome-$lineage"
  label "long_job"

  // conda "$baseDir/conda-envs/busco-env.yaml"

  publishDir "${params.outdir}/busco", mode: 'copy'

  input:
    path genome
    val lineage
    val mode

  output:
    path "busco_${genome.simpleName}_${lineage}/short_summary*", emit: summary

  script:
    def busco_long = params.busco_long ? "--long" : ""

    """
    busco \
    --in $genome \
    --lineage_dataset $lineage \
    --out busco_${genome.simpleName}_${lineage} \
    --mode $mode \
    $busco_long \
    --cpu ${task.cpus}
    """
}
