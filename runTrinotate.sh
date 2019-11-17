
cd /global/scratch/cmt2/transcriptome_raw_data/transcriptomes_Bomarea/bom_mult/trinity
/global/home/users/cmt2/trinityrnaseq-Trinity-v2.8.4/util/support_scripts/get_Trinity_gene_to_trans_map.pl Trinity.fasta > Trinity.fasta.gene_trans_map


module load blast
module load hmmer 
module load sqlite

/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/auto/autoTrinotate.pl  --Trinotate_sqlite /global/scratch/cmt2/transcriptome_raw_data/transcriptomes_Jesus/allium/trinity/Trinotate.sqlite --transcripts Trinity.fasta --gene_to_trans_map Trinity.fasta.gene_trans_map --conf /global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/auto/conf.txt --CPU 10
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite init --gene_trans_map Trinity.fasta.gene_trans_map --transcript_fasta Trinity.fasta --transdecoder_pep Trinity.fasta.transdecoder.pep
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastx swissprot.blastx.outfmt6 
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_swissprot_blastp swissprot.blastp.outfmt6
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_pfam TrinotatePFAM.out
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_rnammer Trinity.fasta.rnammer.gff
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_tmhmm tmhmm.out
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite LOAD_signalp signalp.out
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/Trinotate Trinotate.sqlite report > Trinotate.xls
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/util/extract_GO_assignments_from_Trinotate_xls.pl --Trinotate_xls Trinotate.xls -T -I > Trinotate.xls.gene_ontology
/global/home/users/cmt2/Trinotate-Trinotate-v3.1.1/util/annotation_importer/import_transcript_names.pl Trinotate.sqlite Trinotate.xls
