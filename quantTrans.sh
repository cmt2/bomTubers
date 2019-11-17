cd /global/scratch/cmt2/trans_diff_analysis/
mkdir output
module load bowtie2/2.3.4.1
module load samtools/1.8

trinityrnaseq-Trinity-v2.8.4/util/align_and_estimate_abundance.pl --transcripts Trinity.fasta --est_method RSEM --aln_method bowtie2 --output_dir /output --prep_reference --trinity_mode --samples_file replicate_information.txt --seqType fq  
