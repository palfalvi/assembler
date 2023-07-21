process freebayes_call {
  tag "$contig_index"
  label 'large_mem'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  //publishDir "${params.outdir}/freebayes_polish", mode: 'copy'

  input:
    path contig
    val avg_depth
    path bam
    path baidx
    val contig_index

  output:
    path "*.bcf"

  script:
    """
    samtools faidx ${contig}

    coverage=`printf "%.0f" ${avg_depth}`

    LEN=`wc -l ${contig}.fai | awk '{print \$1}'`

    for j in \$(seq ${contig_index} 100 \$LEN )
    do
      contig=`sed -n \${j}p ${contig}.fai | awk '{print \$1}'`
      contig_no_pipe=`echo \$contig | sed 's/|/_/g'`
      end=`sed -n \${j}p ${contig}.fai | awk '{print \$2}'`

      freebayes --bam $bam --region=\$contig:1-\$end --skip-coverage \$((\$coverage*12)) -f ${contig} | bcftools view --no-version -Ob > \${contig_no_pipe}.bcf
    done
    """
}
