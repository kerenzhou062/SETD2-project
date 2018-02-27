#!/bin/sh

### peak fold_enrichment
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/venn
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/venn
cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/venn
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($14<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' ../xls/allPeak/HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shSetD2_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($15<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' ../xls/allPeak/HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shM14_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($16<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' ../xls/allPeak/HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shM3_FC_gene.txt
awk 'BEGIN{OFS="\t";};{ if(NR>1){ cutoff=1/2;if($17<cutoff){ gsub(/\.[0-9]+/,"" ,$4);print $4"-"NR; } } }' ../xls/allPeak/HepG2_m6a_peak_FC.txt | sort | uniq > ./HepG2_shWTAP_FC_gene.txt
