include { minimap2 } from '../modules/minimap.nf'
include { minimap2_genomes } from '../modules/minimap_genomes.nf'
include { pd_cov } from '../modules/pd_cov.nf'
include { purge_dups } from '../modules/purge_dups.nf'

workflow vgp_polish_wf {
    take:
    genome
    path long_reads
    val platform // map-ont (ONT) map-pb (PB CLR) or asm20 (HiFi)

    main:

    // Coverage calculation
    minimap2( long_reads, genome )
    pd_cov( minimap2.out.map )

    // minimap genome vs genome
    minimap2_genomes( genome, genome )

    // purge_dups
    purge_dups( pd_cov.out.base_cov, pd_cov.out.cutoffs, minimap2_genomes.out.map, genome )

    emit:
        assembly = purge_dups.out.assembly
}
