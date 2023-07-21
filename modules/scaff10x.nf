process scaffX {

  label "assembly"
  publishDir "${params.outdir}/scaff10x", mode: 'copy'

  input:
    path assembly
    tuple val(sample_id), file(reads)

  output:
    path '*.scaff10x.fa', emit: assembly

  script:
    """
    cp ${reads[0]} ${reads[0].simpleName}_data1.fastq.gz
    cp ${reads[1]} ${reads[1].simpleName}_data2.fastq.gz

    touch reads.dat && \
    echo q1=${reads[0].simpleName}_data1.fastq.gz > reads.dat && \
    echo q2=${reads[1].simpleName}_data2.fastq.gz >> reads.dat

    scaff10x \
        -nodes ${task.cpus} \
        -longread 1 \
        -gap 100 \
        -matrix 2000 \
        -reads 10 \
        -score 20 \
        -link 8 \
       -block 50000 \
       -data reads.dat \
        $assembly \
        ${assembly.simpleName}.scaff10x.fa
    """
}
