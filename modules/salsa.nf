process salsa {
  label "small_job"

  conda "$baseDir/conda-envs/salsa-env.yaml"

  publishDir "${params.outdir}/HiC", mode: 'copy'

  input:
    path assembly
    path bam
    path baidx

  output:
    path 'scaffolds/scaffolds_FINAL.fasta', emit: assembly

  script:
    """
      bamToBed -i $bam > ${bam.simpleName}.bed
      sort -k 4 ${bam.simpleName}.bed > tmp && mv tmp ${bam.simpleName}.bed

      samtools faidx $assembly

      run_pipeline.py \
        -a ${assembly} \
        -l ${assembly}.fai \
        -b ${bam.simpleName}.bed \
        -e GATC,GANTC \
        -o scaffolds \
        -m yes
    """
}
