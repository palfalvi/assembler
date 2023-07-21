include { racon } from '../modules/racon.nf'
include { minimap2} from '../modules/minimap.nf'

workflow racon_wf {
    take:
      assembly
      fastq

    main:
    minimap2(fastq, assembly)
    racon(fastq, minimap2.out.map, assembly)



    emit:
        assembly = racon.out.assembly
}
