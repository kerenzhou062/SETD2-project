#!/bin/sh
########HepG2

./HepG2_m6a_seq_FC.sh
./HepG2_m6a_seq_FC_cumulative.sh
./HepG2_m6a_seq_bin_FC.sh
./HepG2_m6a_seq_FC_Venn.sh
./HepG2_m6a_seq_genePeak.sh
./HepG2_m6a_seq_boxplot.sh
./HepG2_m6a_seq_scatter.sh
./HepG2_m6a_seq_exomePeakDE.sh
./HepG2_m6a_seq_motif.sh
./HepG2_m6a_seq_circos_bin.sh
#./HepG2_m6a_seq_tdf.sh
#./HepG2_m6a_seq_stats.sh
./HepG2_m6a_seq_heterochromatin.sh


########Hela

./Hela_m6a_seq_FC.sh
./Hela_m6a_seq_FC_cumulative.sh
./Hela_m6a_seq_bin_FC.sh
./Hela_m6a_seq_FC_Venn.sh
./Hela_m6a_seq_genePeak.sh
./Hela_m6a_seq_boxplot.sh
./Hela_m6a_seq_scatter.sh
./Hela_m6a_seq_motif.sh
./Hela_m6a_seq_circos_bin.sh
#./Hela_m6a_seq_tdf.sh
#./Hela_m6a_seq_stats.sh

########mESCs

./mESCs_m6a_seq_FC.sh
./mMEFs_m6a_seq_bin_FC.sh
./mESCs_m6a_seq_motif.sh
./mESCs_m6a_seq_circos_bin.sh
#./mESCs_m6a_seq_tdf.sh

########mMEFs
./mMEFs_m6a_seq_FC.sh
./mMEFs_m6a_seq_bin_FC.sh
./mMEFs_m6a_seq_motif.sh
#./mMEFs_m6a_seq_tdf.sh
#./mESCs_m6a_seq_circos_bin.sh

