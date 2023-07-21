process hypo {
  label "assembly"

  conda "$baseDir/conda-envs/hypo-env.yaml"

  publishDir "${params.outdir}/short_polished", mode: 'copy'

  input:
    path genome
    tuple val(sample_id), file(short_reads)
    path long_reads
    val genome_size
    path short_bam
    path short_baidx
    val coverage

  output:
    path "*hypo.fasta", emit: assembly

  script:

    """
    # short read
    # minimap2 --secondary=no --MD -L -t ${task.cpus} -ax sr $genome ${short_reads[0]} ${short_reads[1]} | samtools sort -@ $task.cpus -O BAM - > short.bam
    # samtools index short.bam

    #long (ONT)
    minimap2 --secondary=no --MD -L -t ${task.cpus} -ax map-ont $genome ${long_reads} | samtools sort -@ $task.cpus -O BAM - > long.bam
    samtools index long.bam

    # avg_depth=`samtools depth short.bam  |  awk '{sum+=\$3} END { print sum/NR}'`

    echo -e "${short_reads[0]}\n${short_reads[1]}" > names.txt

    hypo \
    --draft $genome \
    --reads-short @names.txt \
    --size-ref $genome_size \
    --coverage-short $coverage \
    --processing-size 96 \
    --bam-sr $short_bam \
    --bam-lr long.bam \
    --threads ${task.cpus} \
    --output ${genome.simpleName}_hypo.fasta
    """
}
