include { minimap2 } from '../modules/minimap.nf'
include { minimap2_genomes } from '../modules/minimap_genomes.nf'
include { pd_cov } from '../modules/pd_cov.nf'
include { purge_dups } from '../modules/purge_dups.nf'

workflow purge_dups_wf {
    take:
    assembly
    long_reads

    main:

    // Coverage calculation
    minimap2( long_reads, assembly )
    pd_cov( minimap2.out.map )

    // minimap genome vs genome
    minimap2_genomes( assembly, assembly )

    // purge_dups
    purge_dups( pd_cov.out.base_cov, pd_cov.out.cutoffs, minimap2_genomes.out.map, assembly )

    emit:
        assembly = purge_dups.out.assembly
}
