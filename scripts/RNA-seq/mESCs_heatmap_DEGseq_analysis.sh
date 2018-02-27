#!/bin/sh

### mESCs, RNA-seq
cd /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mESCs/heatmap/
annotation=/data/zhoukr/hhl_setd2_m6a/reference/annotation/gencode.vM9.annotation.all.longest.exon.bed6
DEGseqFile=/data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/mouse/DEGseq/mouse_RNA-seq_DEGseq_intersected_FC.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==3{
  geneName = hashGeneIdName[$1]
  EXP[geneName]= $2"\t"$3"\t"$4"\t"$5;
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (i in EXP) {
      print i, EXP[i];
    }
  }
}
' $annotation mmu04550_stem.txt $DEGseqFile > heatmapData/DEGseq/mESCs_DEGseq_mmu04550_stem_heatmap_log2FC.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==3{
  geneName = hashGeneIdName[$1]
  EXP[geneName]= $2"\t"$3"\t"$4"\t"$5;
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (i in EXP) {
      print i, EXP[i];
    }
  }
}
' $annotation mmu04550_differ.txt $DEGseqFile > heatmapData/DEGseq/mESCs_DEGseq_mmu04550_differ_heatmap_log2FC.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==3{
  geneName = hashGeneIdName[$1]
  EXP[geneName]= $2"\t"$3"\t"$4"\t"$5;
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (i in EXP) {
      print i, EXP[i];
    }
  }
}
' $annotation CONRAD_STEM_CELL_18849962.txt $DEGseqFile > heatmapData/DEGseq/mESCs_DEGseq_18849962_heatmap_log2FC.txt

awk 'BEGIN{OFS="\t";}
ARGIND==1{
  split($4,geneIdI,":");
  geneIdi=geneIdI[1];
  geneNamei=geneIdI[2];
  hashGeneIdName[geneIdi]=geneNamei;
}
ARGIND==2{
  Ulower = "";
  if(length($1)==1){
    Ulower=toupper($1);
  }else{
    Ulower=toupper(substr($1,1,1)) "" tolower(substr($1,2));
  }
  hashGene[Ulower] = 1;
}
ARGIND==3{
  geneName = hashGeneIdName[$1]
  EXP[geneName]= $2"\t"$3"\t"$4"\t"$5;
}
END{
  print "GeneName\tshContD0\tshContD6\tshSETD2D0\tshSetD2D6";
  for (i in hashGene) {
    if (i in EXP) {
      print i, EXP[i];
    }
  }
}
' $annotation MIKKELSEN_PLURIPOTENT_STATE_UP_18509334.txt $DEGseqFile > heatmapData/DEGseq/mESCs_DEGseq_18509334_heatmap_log2FC.txt
