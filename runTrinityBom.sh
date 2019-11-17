cd /global/scratch/cmt2/transcriptome_raw_data/transcriptomes_Bomarea/bom_mult

Trinity --SS_lib_type RF --seqType fq --max_memory 400G --left bom_mult_1.fq --right bom_mult_2.fq --CPU 20 --output trinity

exit