mv CHe-TW-HL-A_S1_L008_R1_001.fastq.gz Hela_m6A-seq_shCont_input_rep1.fastq.gz
mv CHe-TW-HL-C_S2_L008_R1_001.fastq.gz Hela_m6A-seq_shM3_input_rep1.fastq.gz
mv CHe-TW-HL-D_S3_L008_R1_001.fastq.gz Hela_m6A-seq_shM3_input_rep2.fastq.gz
mv CHe-TW-HL-E_S4_L008_R1_001.fastq.gz Hela_m6A-seq_shM14_input_rep1.fastq.gz
mv CHe-TW-HL-F_S5_L008_R1_001.fastq.gz Hela_m6A-seq_shM14_input_rep2.fastq.gz
mv CHe-TW-HL-G_S6_L008_R1_001.fastq.gz Hela_m6A-seq_shWTAP_input_rep1.fastq.gz
mv CHe-TW-HL-K_S7_L008_R1_001.fastq.gz Hela_m6A-seq_shCont_IP_rep1.fastq.gz
mv CHe-TW-HL-M_S8_L008_R1_001.fastq.gz Hela_m6A-seq_shM3_IP_rep1.fastq.gz
mv CHe-TW-HL-N_S9_L008_R1_001.fastq.gz Hela_m6A-seq_shM3_IP_rep2.fastq.gz
mv CHe-TW-HL-O_S10_L008_R1_001.fastq.gz Hela_m6A-seq_shM14_IP_rep1.fastq.gz
mv CHe-TW-HL-P_S11_L008_R1_001.fastq.gz Hela_m6A-seq_shM14_IP_rep2.fastq.gz
mv CHe-TW-HL-Q_S12_L008_R1_001.fastq.gz Hela_m6A-seq_shWTAP_IP_rep1.fastq.gz
mv CHe-TWu-14S-TW-HL-B_S1_L008_R1_001.fastq.gz Hela_m6A-seq_shCont_input_rep2.fastq.gz
mv CHe-TWu-14S-TW-HL-L_S2_L008_R1_001.fastq.gz Hela_m6A-seq_shCont_IP_rep2.fastq.gz
mv CHe-TWu-14S-TW-HL-H_S3_L008_R1_001.fastq.gz Hela_m6A-seq_shWTAP_input_rep2.fastq.gz
mv CHe-TWu-14S-TW-HL-I_S4_L008_R1_001.fastq.gz Hela_m6A-seq_shSetD2_input_rep1.fastq.gz
mv CHe-TWu-14S-TW-HL-J_S5_L008_R1_001.fastq.gz Hela_m6A-seq_shSetD2_input_rep2.fastq.gz
mv CHe-TWu-14S-TW-HL-R_S6_L008_R1_001.fastq.gz Hela_m6A-seq_shWTAP_IP_rep2.fastq.gz
mv CHe-TWu-14S-TW-HL-S_S7_L008_R1_001.fastq.gz Hela_m6A-seq_shSetD2_IP_rep1.fastq.gz
mv CHe-TWu-14S-TW-HL-T_S8_L008_R1_001.fastq.gz Hela_m6A-seq_shSetD2_IP_rep2.fastq.gz


gzip -dq Hela_m6A-seq_shCont_input_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM3_input_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM3_input_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shM14_input_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM14_input_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shWTAP_input_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shCont_IP_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM3_IP_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM3_IP_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shM14_IP_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shM14_IP_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shWTAP_IP_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shCont_input_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shCont_IP_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shWTAP_input_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shSetD2_input_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shSetD2_input_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shWTAP_IP_rep2.fastq.gz &
gzip -dq Hela_m6A-seq_shSetD2_IP_rep1.fastq.gz &
gzip -dq Hela_m6A-seq_shSetD2_IP_rep2.fastq.gz &
