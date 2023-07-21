include { bwa_index } from '../modules/bwa_index.nf'
include { bwa_mem } from '../modules/bwa_mem.nf'
include { bam_coverage } from '../modules/bam_coverage.nf'
include { bam_merge } from '../modules/bam_merge.nf'
include { hypo } from '../modules/hypo.nf'

workflow hypo_wf {
    take:
      assembly
      short_r
      long_r
      genome_size

    main:
    // Mapping with bwa-mem

    bwa_index( assembly, "" )

    bwa_mem( short_r, assembly, bwa_index.out.index )

    bam_merge( bwa_mem.out.bam.collect(), "short_reads" )

    bam_coverage( bam_merge.out.bam )

    bam_coverage.out.coverage.subscribe { println "Short read coverage is ${it}x." }

    hypo( assembly, short_r, long_r, genome_size, bam_merge.out.bam, bam_merge.out.baidx, bam_coverage.out.coverage )



    emit:
        assembly = hypo.out.assembly
}
