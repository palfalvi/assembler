process purge_dups {

  label 'long_job'

  conda "$baseDir/conda-envs/purge_dups-env.yaml"

  publishDir "${params.outdir}/purge_dups", mode: 'copy'

  input:
  path base_cov
  path cutoffs
  path self_alignment
  assembly

  output:
    path '*_purged.fasta', emit: assembly
    path '*_hap.fasta',    emit: haplo

  script:
    """
    purge_dups -2 -T $cutoffs $base_cov ${genome.simpleName}_self.paf.gz > dups.bed 2> purge_dups.log

    get_seqs dups.bed $assembly > ${genome.simpleName}_purged.fasta 2> ${genome.simpleName}_hap.fasta
    """
}
