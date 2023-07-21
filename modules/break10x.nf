process breakX {
  tag "$assembly"
  label "small_plus"
  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.break10x.fa', emit: assembly

  script:
    """
    SCAFF10X=`which break10x`

    cp ${reads[0]} ${reads[0].simpleName}_data1.fastq.gz
    cp ${reads[1]} ${reads[1].simpleName}_data2.fastq.gz

    touch reads.dat && \
    echo q1=${reads[0].simpleName}_data1.fastq.gz > reads.dat && \
    echo q2=${reads[1].simpleName}_data2.fastq.gz >> reads.dat

    \$SCAFF10X \
        -nodes ${task.cpus} \
        -data reads.dat \
        $assembly \
        ${assembly.simpleName}.break10x.fa \
        ${assembly.simpleName}.break10x.name
    """
}
