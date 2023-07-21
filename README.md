# Genome Assembly Pipeline

Test on Drosera rotundifolia polishing:

nohup \
nextflow run /mnt/gpfsA/home/peegee/nf/assembler -profile cde -bg \
--assembly '/mnt/gpfsB/scratch/peegee/drot_gala/assemblies/*.fasta' \
--genome_size 1.15g \
--purge \
--polish_r10 \
--fastq /mnt/gpfsB/scratch/peegee/drot_necat/ont_reads/raw/Drot_guppy5.all.fastq \
--fastq_r9 /mnt/gpfsB/scratch/peegee/drot_necat/ont_reads/raw/Drot_R9_guppy.all.fastq \
--fastq_r10 /mnt/gpfsB/scratch/peegee/drot_necat/ont_reads/raw/Drot_R10_guppy.all.fastq \
--short_polish hypo \
--short_reads '/mnt/gpfsB/scratch/peegee/drot_necat/drot_necat2/8-scaff10x/drot_10XG-BC_{1,2}.fastq.gz' \
&
