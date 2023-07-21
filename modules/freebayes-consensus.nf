process freebayes_consensus {

  label 'assembly'

  conda "$baseDir/conda-envs/freebayes-env.yaml"

  publishDir "${params.outdir}/short_polished", mode: 'copy'

  input:
    path assembly
    path bcf_list

  output:
    path "*freebayes.fasta", emit: assembly

  script:
    """

    for file in ./tarseq_*.bcf; do bcftools index \$file  -f; done

    ls *.bcf > bcf_files.txt

    bcftools concat -af bcf_files.txt | bcftools view -Ou -e'type="ref"' --threads ${task.cpus} | bcftools norm --threads ${task.cpus} -Ob -f ${assembly} -o ${assembly.simpleName}.bcf
    bcftools index ${assembly.simpleName}.bcf

    bcftools consensus -i'QUAL>1 && (GT="AA" || GT="Aa")' -Hla -f ${assembly} ${assembly.simpleName}.bcf > ${assembly.simpleName}_freebayes.fasta
    """
}
