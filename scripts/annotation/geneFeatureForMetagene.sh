#!/bin/sh
### gene body
#awk '
#BEGIN{OFS="\t";FS="\t";exLen=10000;}
#ARGIND==1{
#    chr=$1;size=$2;
#    hash[chr]=size;
#}
#ARGIND==2{
#    if(FNR>5){
#       if($3=="gene"){
#           if($9~/gene_type "protein_coding"/){
#               chr=$1;start=$4-1;end=$5;strand=$7;
#               match($9,/gene_id "(\w+\.?\w+)";/,pattern);geneID=pattern[1];
#               match($9,/gene_name "(\w+\.?\w+)";/,pattern);geneName=pattern[1];
#               P5Start=start-exLen;if(P5Start<0){P5Start=0;};
#               if(strand=="+"){name=geneID":"geneName"|5PStart|1";}else{name=geneID":"geneName"|3PEnd|1";}
#               print chr,P5Start,start,name,0,strand;
#               name=geneID":"geneName"|geneBody|1";
#               print chr,start,end,name,0,strand;
#               P3End=end+exLen;if(P3End>hash[chr]){P3End=hash[chr]};
#               if(strand=="+"){name=geneID":"geneName"|3PEnd|1";}else{name=geneID":"geneName"|5PStart|1";}
#               print chr,end,P3End,name,0,strand;
#           }
#       }
#   }
#}' human_hg19_chrsize.txt gencode.v24lift37.annotation.gtf \
#   | sort -t $'\t' -k1,1V -k2,2n > gencode.v24lift37.mRNA.geneFeature.bed6
#
#awk '
#BEGIN{OFS="\t";FS="\t";exLen=10000;}
#ARGIND==1{
#    chr=$1;size=$2;
#    hash[chr]=size;
#}
#ARGIND==2{
#    if(FNR>5){
#       if($3=="gene"){
#           if($9~/gene_type "protein_coding"/){
#               chr=$1;start=$4-1;end=$5;strand=$7;
#               match($9,/gene_id "(\w+\.?\w+)";/,pattern);geneID=pattern[1];
#               match($9,/gene_name "(\w+\.?\w+)";/,pattern);geneName=pattern[1];
#               P5Start=start-exLen;if(P5Start<0){P5Start=0;};
#               if(strand=="+"){name=geneID":"geneName"|5PStart|1";}else{name=geneID":"geneName"|3PEnd|1";}
#               print chr,P5Start,start,name,0,strand;
#               name=geneID":"geneName"|geneBody|1";
#               print chr,start,end,name,0,strand;
#               P3End=end+exLen;if(P3End>hash[chr]){P3End=hash[chr]};
#               if(strand=="+"){name=geneID":"geneName"|3PEnd|1";}else{name=geneID":"geneName"|5PStart|1";}
#               print chr,end,P3End,name,0,strand;
#           }
#       }
#   }
#}' mouse_mm10_chrsize.txt gencode.vM9.annotation.gtf \
#   | sort -t $'\t' -k1,1V -k2,2n > gencode.vM9.mRNA.geneFeature.bed6

### cds start and cds end
awk '
BEGIN{OFS="\t";FS="\t";exLen=5000;}
ARGIND==1{
    chr=$1;size=$2;
    hash[chr]=size;
}
ARGIND==2{
    if(FNR>5){
      chr=$1;split($4,nameArr,"|");gene=nameArr[1];strand=$6;thickStart=$7;thickEnd=$8;
      if(strand=="+"){
        P5Start=thickStart-exLen;if(P5Start<0){P5Start=0;};
        P3End=thickEnd+exLen;if(P3End>hash[chr]){P3End=hash[chr]};
        name=gene"|5PStart|1";
        print chr,P5Start,thickStart,name,0,strand;
        name=gene"|geneBody|1";
        print chr,thickStart,thickEnd,name,0,strand;
        name=gene"|3PEnd|1";
        print chr,thickEnd,P3End,name,0,strand;
      }else{
        P5Start=thickEnd+exLen;if(P5Start>hash[chr]){P5Start=hash[chr]};
        P3End=thickStart-exLen;if(P3End<0){P3End=0;};
        name=gene"|5PStart|1";
        print chr,thickEnd,P5Start,name,0,strand;
        name=gene"|geneBody|1";
        print chr,thickStart,thickEnd,name,0,strand;
        name=gene"|3PEnd|1";
        print chr,P3End,thickStart,name,0,strand;
      }
    }
}' human_hg19_chrsize.txt gencode.v24lift37.annotation.mRNA.longest.exon.bed12 \
   | sort -t $'\t' -k1,1V -k2,2n > gencode.v24lift37.mRNA.geneFeature.bed6

awk '
BEGIN{OFS="\t";FS="\t";exLen=5000;}
ARGIND==1{
    chr=$1;size=$2;
    hash[chr]=size;
}
ARGIND==2{
    if(FNR>5){
      chr=$1;split($4,nameArr,"|");gene=nameArr[1];strand=$6;thickStart=$7;thickEnd=$8;
      if(strand=="+"){
        P5Start=thickStart-exLen;if(P5Start<0){P5Start=0;};
        P3End=thickEnd+exLen;if(P3End>hash[chr]){P3End=hash[chr]};
        name=gene"|5PStart|1";
        print chr,P5Start,thickStart,name,0,strand;
        name=gene"|geneBody|1";
        print chr,thickStart,thickEnd,name,0,strand;
        name=gene"|3PEnd|1";
        print chr,thickEnd,P3End,name,0,strand;
      }else{
        P5Start=thickEnd+exLen;if(P5Start>hash[chr]){P5Start=hash[chr]};
        P3End=thickStart-exLen;if(P3End<0){P3End=0;};
        name=gene"|5PStart|1";
        print chr,thickEnd,P5Start,name,0,strand;
        name=gene"|geneBody|1";
        print chr,thickStart,thickEnd,name,0,strand;
        name=gene"|3PEnd|1";
        print chr,P3End,thickStart,name,0,strand;
      }
    }
}' mouse_mm10_chrsize.txt gencode.vM9.annotation.mRNA.longest.exon.bed12 \
   | sort -t $'\t' -k1,1V -k2,2n > gencode.vM9.mRNA.geneFeature.bed6
