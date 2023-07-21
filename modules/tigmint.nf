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
    tigmint-make arcs draft=myassembly reads=myreads
    """
}
