include { bwa_index as bwa_idx_10x } from '../modules/bwa_index.nf'
include { bwa_mem_hic as bwa_mem_hic1; bwa_mem_hic as bwa_mem_hic2 } from '../modules/bwa_mem_hic.nf'
//include { bwa_mem as bwa_mem_hic2 } from '../modules/bwa_mem.nf'
include { arima_filter as filter_5ends1; arima_filter as filter_5ends2 } from '../modules/arima_filter.nf'
//include { arima_filter as filter_5ends2 } from '../modules/arima_filter.nf'
include { arima_qc } from '../modules/arima_qc.nf'
include { arima_add_read_group as add_read_group } from '../modules/arima_add_read_group.nf'
include { picard_mark_duplicates as mark_duplicates } from '../modules/picard_mark_duplicates.nf'
include { arima_stats as calc_stats } from '../modules/arima_stats.nf'
include { samtools_index} from '../modules/samtools_index.nf'


workflow arima_mapping {
    take:
      genome
      reads

    main:
      bwa_idx_10x(genome, "-a bwtsw")

      bwa_mem_hic1( reads, 0, genome, bwa_idx_10x.out.index )
      bwa_mem_hic2( reads, 1, genome, bwa_idx_10x.out.index )

      filter_5ends1(bwa_mem_hic1.out.bam)
      filter_5ends2(bwa_mem_hic2.out.bam)

      arima_qc(filter_5ends1.out.bam, filter_5ends2.out.bam, genome)

      add_read_group(arima_qc.out.bam)

      mark_duplicates(add_read_group.out.bam)

      samtools_index(mark_duplicates.out.bam, "")

      calc_stats(mark_duplicates.out.bam, samtools_index.out.baidx)


    emit:
        bam = mark_duplicates.out.bam
        baidx = samtools_index.out.baidx
        stats = calc_stats.out.stats
}
