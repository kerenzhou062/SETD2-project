#!/bin/sh
########HepG2

./HepG2_chip-seq_FC.sh
./HepG2_chip-seq_bin.sh
./HepG2_chip-seq_venn.sh
./HepG2_featureCount.sh
./HepG2_chip-seq_genePeak.sh
./HepG2_chip-seq_circos.sh

########Hela
./Hela_chip-seq_FC.sh
./Hela_chip-seq_bin.sh
./Hela_chip-seq_venn.sh
./Hela_chip-seq_genePeak.sh
./Hela_featureCount.sh
./Hela_chip-seq_circos.sh

########mESCs
./mESCs_chip-seq_bin.sh

########mMEFs
./mMEFs_chip-seq_bin.sh
