#!/bin/sh

### mESCs, log2(FPKM+0.01)
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/heatmap/
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/RSEM

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, log(EXPa[i]+0.01)/log(2), log(EXPc[i]+0.01)/log(2), log(EXPb[i]+0.01)/log(2), log(EXPd[i]+0.01)/log(2);
    }
  }
}
' mmu04550_stem.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_mmu04550_stem_heatmap_log2.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, log(EXPa[i]+0.01)/log(2), log(EXPc[i]+0.01)/log(2), log(EXPb[i]+0.01)/log(2), log(EXPd[i]+0.01)/log(2);
    }
  }
}
' mmu04550_differ.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_mmu04550_differ_heatmap_log2.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, log(EXPa[i]+0.01)/log(2), log(EXPc[i]+0.01)/log(2), log(EXPb[i]+0.01)/log(2), log(EXPd[i]+0.01)/log(2);
    }
  }
}
' CONRAD_STEM_CELL_18849962.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_18849962_heatmap_log2.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, log(EXPa[i]+0.01)/log(2), log(EXPc[i]+0.01)/log(2), log(EXPb[i]+0.01)/log(2), log(EXPd[i]+0.01)/log(2);
    }
  }
}
' MIKKELSEN_PLURIPOTENT_STATE_UP_18509334.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_18509334_heatmap_log2.txt

### mESCs, FPKM
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/heatmap/
geneExpDir=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/RSEM

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, EXPa[i], EXPc[i], EXPb[i], EXPd[i];
    }
  }
}
' mmu04550_stem.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_mmu04550_stem_heatmap.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, EXPa[i], EXPc[i], EXPb[i], EXPd[i];
    }
  }
}
' mmu04550_differ.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_mmu04550_differ_heatmap.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, EXPa[i], EXPc[i], EXPb[i], EXPd[i];
    }
  }
}
' CONRAD_STEM_CELL_18849962.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_18849962_heatmap.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==2{
  if(FNR>1){
    EXP[$2] = 1;
    EXPa[$2]= $3;
    EXPb[$2]= $4;
  }
}
ARGIND==3{
  if(FNR>1){
    EXP[$2] += 1;
    EXPc[$2]= $3;
    EXPd[$2]= $4;
  }
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (EXP[i] == 2){
      print i, EXPa[i], EXPc[i], EXPb[i], EXPd[i];
    }
  }
}
' MIKKELSEN_PLURIPOTENT_STATE_UP_18509334.txt $geneExpDir/mESCs_RNA-seq_RSEM_D0_input.txt $geneExpDir/mESCs_RNA-seq_RSEM_D6_input.txt > heatmapData/fpkm/mESCs_RSEM_input_18509334_heatmap.txt
