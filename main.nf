#!/usr/bin/env nextflow
nextflow.enable.dsl=2

def helpMessage() {
log.info """
=======================================================
           __
         /    )
--------/---------__----__----__---_--_----__----------
       /  --,   /___) /   ) /   ) / /  ) /___)
______(____/___(___ _/___/_(___/_/_/__/_(___ __________
    __
    / |                           /     /
---/__|---__---__----__---_--_---/__---/----__---)__-
  /   |  (_ ` (_ ` /___) / /  ) /   ) /   /___) /   )
_/____|_(__)_(__)_(___ _/_/__/_(___/_/___(___ _/_____

=======================================================

Usage:


Mandatory arguments:
      --fastq                        Long read fastq file
      --genome_size                  Expected size of genome.

Universal arguments
      --outdir                       Output directory name. [results]
      -profile                       Sets the running environment. Default is NIBB-BIAS5 PBSPro. 'cde' and 'local' are available to run on NIBB-CDE server or on local machine.
      -bg                            Run the pipeline in the background.
      -resume                        Resume interrupted run and try to catch previously finished processes.



Assembly mode:
      --short                        Short read fastq file. Used for optional polishing.
      --polish                       True/False or software?
      --scaffX                       Scaff10x executable have to be in PATH!


""".stripIndent()
}

params.help = false
if (params.help){
helpMessage()
exit 0
}


// Include assemblers
include { canu } from './modules/canu.nf'
//include { necat } from './modules/necat.nf'
include { raven } from './modules/raven.nf'
include { flye } from './modules/flye.nf'
include { miniasm } from './modules/miniasm.nf'
include { wtdbg } from './modules/wtdbg.nf'

// Include scaffolders
include { scaffX } from './modules/scaff10x.nf'
include { breakX } from './modules/break10x.nf'
include { purge_dups } from './modules/purge_dups.nf'
include { debarcodeX } from './modules/debarcode10x.nf'
include { arima_mapping } from './modules/arima_mapping.nf'
include { salsa } from './modules/salsa.nf'

// Include polishing tools
include { minimap2 as minimap2_1; minimap2 as minimap2_2; minimap2 as minimap2_3 } from './modules/minimap.nf'
include { bwa_index } from './modules/bwa_index.nf'
include { bwa_mem } from './modules/bwa_mem.nf'
include { bam_coverage } from './modules/bam_coverage.nf'
include { bam_merge } from './modules/bam_merge.nf'
include { racon_wf } from './modules/racon_wf.nf'
include { medaka } from './modules/medaka.nf'
include { medaka as medaka_r10 } from './modules/medaka.nf'
include { hypo_wf } from './modules/hypo_wf.nf'
include { vgp_polish_wf } from './modules/vgp_polish_wf.nf'
include { pilon } from './modules/pilon.nf'

// Include QC tools
include { busco } from './modules/busco.nf'
include { quast } from './modules/quast.nf'
include { multiqc } from './modules/multiqc.nf'
include { kat } from './modules/kat.nf'


workflow {

/////////////// ASSEMBLY PIPELINE ///////////////
log.info """
=======================================================
           __
         /    )
--------/---------__----__----__---_--_----__----------
       /  --,   /___) /   ) /   ) / /  ) /___)
______(____/___(___ _/___/_(___/_/_/__/_(___ __________
    __
    / |                           /     /
---/__|---__---__----__---_--_---/__---/----__---)__-
  /   |  (_ ` (_ ` /___) / /  ) /   ) /   /___) /   )
_/____|_(__)_(__)_(___ _/_/__/_(___/_/___(___ _/_____

=======================================================
""".stripIndent()


log.info ">>> Starting Genome Assembler pipeline ... "

  // Canu, out.assembly, out.gfa
  if ( params.canu ) {
    log.info ">>> Starting canu assembly."
  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for canu. Please provide with --genome_size."

    canu(params.fastq, params.genome_size)

    canu_assembly = canu.out.assembly
    canu_gfa = canu.out.gfa
  }

  // Flye: out.assembly, out.gfa
  if ( params.flye ) {
    log.info ">>> Starting flye assembly."
  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for flye. Please provide with --genome_size."

    flye(params.fastq, params.genome_size)

    flye_assembly = flye.out.assembly
    flye_gfa = flye.out.gfa
  }

  // Miniasm out.assembly, out.gfa
  if ( params.miniasm ) {
    log.info ">>> Starting miniasm assembly."
  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    miniasm(params.fastq)

    miniasm_assembly = miniasm.out.assembly
    miniasm_gfa = miniasm.out.gfa
  }

  // wtdbg2: out.assembly
  if ( params.wtdbg2 ) {
    log.info ">>> Starting wtdbg2 assembly."
  //  params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."
  //  params.genome_size ? log.info "Estimated genome size: $it" : error "Estimated genome size is missing but needed for wtdbg2. Please provide with --genome_size."

    wtdbg(params.fastq, params.genome_size)

    wtdbg_assembly = wtdbg.out.assembly
  }

  // RAVEN assembler: out.assembly out.gfa
  if ( params.raven ) {
    log.info ">>> Starting raven assembly."
    // params.fastq ? log.info "Fastq file provided: $it" : error "Fastq file is not provided. Please specify with --fastq parameter."

    raven(params.fastq)

    raven_assembly = raven.out.assembly
    raven_gfa = raven.out.gfa

  }

  // NECAT assembly
  if ( params.necat ) {
    log.info ">>> Starting NECAT assembly."

    //necat(params.fastq)
    //necat_assembly = necat.out.assembly
    //necat_gfa = necat.out.gfa

  }

  if ( params.assembly ) {
    log.info ">>> Primary assembly provided."

    assembly = Channel.fromPath( params.assembly )

  }

//////////////////////////////////
//////// HAPLOTIG PURGING ////////
/////////////////////////////////

  if ( params.purge ) {
    log.info ">>> Purging haplotigs."
    purge_dups(assembly, params.fastq, "map-ont")
    assembly = purge_dups.out.purged
    // purge_duplicates https://depot.galaxyproject.org/singularity/purge_dups:1.2.5--h5bf99c6_1
  }

  /////////////////////////////////
  //////// 10X scaffolding ////////
  /////////////////////////////////

  if ( params.scaffX ) {
    // use scaff10x and break10x
    // need to install https://github.com/wtsi-hpag/Scaff10X manually

    if ( !params.linked_reads) {
      error 'Short reads are not provided. Please provide short reads as --short_reads /path/to/short.fastq or as --linked_reads /path/to/linked.fastq'
    } else {
      linked_r = Channel.fromFilePairs( params.linked_reads )
      linked_r.subscribe {  println "Linked reads provided: $it"  }

      log.info ">>> Scaffolding primary assembly with Scaff10x."

      //scaffX( assembly, linked_r )
      //breakX( scaffX.out.assembly, linked_r )

      // Tigmint + ARCS

      assembly = breakX.out.assembly
  }

      }


 /////////////////////////////////
 //////// HiC scaffolding ////////
 /////////////////////////////////

        if ( params.hic_reads ) {
          // use HiC reads for scaffolding

          hic_r = Channel.fromFilePairs( params.hic_reads )

          hic_r.subscribe {  println "HiC reads provided: $it"  }
          log.info ">>> Scaffolding primary assembly with Salsa."

          //arima_mapping( assembly, hic_r )
          //salsa( assembly, arima_mapping.out.bam, arima_mapping.out.baidx )
          assembly = salsa.out.assembly
        }

/////////////////////////////////////
//////// LONG READ POLISHING ////////
/////////////////////////////////////


  if ( params.polish ) {
    // racon and medaka polishing
    log.info ">>> Polishing assembly with long reads."

    racon_wf(assembly, params.fastq)

    medaka(params.fastq, racon_wf.out.assembly, params.r9_medaka)

    //medaka.out.assembly

    assembly = medaka.out.assembly

  } else if ( params.polish_r10 ) {
      // racon and medaka polishing
      log.info ">>> Polishing assembly with long reads."

      racon_wf(assembly, params.fastq)

      medaka(params.fastq_r9, racon_wf.out.assembly, params.r9_medaka)

      medaka_r10(params.fastq_r10, medaka.out.assembly, params.r10_medaka)

      //medaka.out.assembly

      assembly = medaka_r10.out.assembly

    }else if (params.racon_polish) {
    // only racon polishing
    log.info ">>> Polishing assembly with long reads."

    racon_wf(params.fastq, assembly)

    assembly = racon_wf.out.assembly
  } else if (params.medaka_polish) {
    // only medaka polishing
    log.info ">>> Polishing assembly with long reads."

    medaka(params.fastq, assembly)

    assembly = medaka.out.assembly
  }

//////////////////////////////////////
//////// SHORT READ POLISHING ////////
//////////////////////////////////////

  if ( params.short_polish ) {

    if ( !params.short_reads && !params.linked_reads) {
      error 'Short reads are not provided. Please provide short reads as --short_reads /path/to/short.fastq or as --linked_reads /path/to/linked.fastq'
    } else {
      // Short read mapping
      if ( params.short_reads ) {

        short_r = Channel.fromFilePairs( params.short_reads )
        short_r.subscribe {  println "Short reads provided: $it"  }
        log.info ">>> Polishing assembly with short reads."

      } else if ( params.linked_reads ){

        short_r = Channel.fromFilePairs( params.linked_reads )
        short_r.subscribe {  println "Linked reads provided: $it"  }
        log.info ">>> Polishing assembly with linked reads."

        // Use scaff10x to debarcode reads

        }

      }

    if ( params.short_polish == 'freebayes' | params.short_polish == 'vgp' | params.short_polish == true ) {
      // Verterae Genome Project polishing pipeline with freebayes and bcftools

      vgp_polish_wf( assembly, short_r )

      assembly = vgp_polish_wf.out.assembly

    } else if ( params.short_polish == 'nextpolish' ) {
      // Nextpolish polishing
      if ( !params.nextpolish ) {
        error 'NextPolish source is not provided. Please install NextPolish locally and provide as --nextpolish /PATH/TO/nextpolish or consider using another polishing method (e.g. hypo).'
      }

      error 'Sorry. Nextpolish is not yet implemented.'

    } else if ( params.short_polish == 'pilon' ) {
      // pilon polishing
      error 'Sorry. PILON is not yet implemented.'
      pilon( assembly, short_bam, short_baidx )

      polished_assembly = pilon.out.assembly

    }  else if ( params.short_polish == 'hypo' ) {
      // HyPo polishing

      hypo_wf( assembly, short_r, params.fastq, params.genome_size )

      assembly = hypo_wf.out.assembly
    }

  }





/////////////////////////////
//////// ASSEMBLY QC ////////
/////////////////////////////


  if ( !params.skip_qc) {
    // QC
    quast(assembly)

    busco(assembly, Channel.fromList(params.busco_lineages), "genome")

    multiqc(quast.out.summary.mix(busco.out).collect(), "$baseDir/${params.outdir}")

    if ( params.short_reads ) {

      short_reads = Channel.fromFilePairs( params.short_reads )

      kat(assembly, short_reads)

    }
  }

}


workflow.onComplete {
    if ( workflow.success ) {
      log.info "[$workflow.complete] >> Assembler pipeline finished SUCCESSFULLY after $workflow.duration ."
      log.info "[$workflow.complete] >> "
      log.info "[$workflow.complete] >> You can find further help on https://github.com/palfalvi"
    } else {
      log.info "[$workflow.complete] >> The script quit with ERROR after ${workflow.duration}."
      log.info "[$workflow.complete] >> Please revise your code and resubmit jobs with the -resume option or reach out for help at https://github.com/palfalvi."
    }
}
