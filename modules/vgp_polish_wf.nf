include { bwa_index } from '../modules/bwa_index.nf'
include { bwa_mem } from '../modules/bwa_mem.nf'
include { bam_coverage } from '../modules/bam_coverage.nf'
include { bam_merge } from '../modules/bam_merge.nf'
include { freebayes_call } from '../modules/freebayes-call.nf'
include { freebayes_consensus } from '../modules/freebayes-consensus.nf'

workflow vgp_polish_wf {
    take:
      assembly
      short_r

    main:
    // Mapping with bwa-mem

    bwa_index( assembly, "" )

    bwa_mem( short_r, assembly, bwa_index.out.index )

    bam_merge( bwa_mem.out.bam.collect(), "short_reads" )

    bam_coverage( bam_merge.out.bam )

    bam_coverage.out.coverage.subscribe { println "Short read coverage is ${it}x." }

    contig_index = Channel.from(1..100) // Split into 100 parallel processes

    freebayes_call( assembly, bam_coverage.out.coverage.first(), bam_merge.out.bam.first(), bam_merge.out.baidx.first(), contig_index )

    // how collectFile works? Concatenate files or colelcts all files?
    //freebayes_call.out.collectFile(name: 'concat_list.txt', newLine: true, sort: true).set { bcf_list }

    freebayes_consensus( assembly, freebayes_call.out.collect() )

    emit:
        assembly = freebayes_consensus.out.assembly
}
