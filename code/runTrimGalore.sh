cd /global/scratch/cmt2/transcriptome_raw_data/transcriptomes_Bomarea
mkdir outputdir36bp

trim_galore -o outputdir36bp --paired CSJM003A_S33_L006_R1_001.fastq.gz CSJM003A_S33_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003B_S34_L006_R1_001.fastq.gz CSJM003B_S34_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003C_S35_L006_R1_001.fastq.gz CSJM003C_S35_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003D_S36_L006_R1_001.fastq.gz CSJM003D_S36_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003E_S37_L006_R1_001.fastq.gz CSJM003E_S37_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003F_S38_L006_R1_001.fastq.gz CSJM003F_S38_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003G_S39_L006_R1_001.fastq.gz CSJM003G_S39_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003H_S40_L006_R1_001.fastq.gz CSJM003H_S40_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003I_S41_L006_R1_001.fastq.gz CSJM003I_S41_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003J_S42_L006_R1_001.fastq.gz CSJM003J_S42_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003L_S44_L006_R1_001.fastq.gz CSJM003L_S44_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
trim_galore -o outputdir36bp --paired CSJM003k_S43_L006_R1_001.fastq.gz CSJM003k_S43_L006_R2_001.fastq.gz --retain_unpaired --dont_gzip --length 36
