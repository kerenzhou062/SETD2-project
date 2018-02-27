#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/circos
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/HepG2_shSetD2_m6a.xls ../xls/HepG2_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > HepG2_m6a_shSetD2_FC_circos.txt
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/HepG2_shM14_m6a.xls ../xls/HepG2_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > HepG2_m6a_shM14_FC_circos.txt
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/HepG2_shM3_m6a.xls ../xls/HepG2_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > HepG2_m6a_shM3_FC_circos.txt
exome2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/HepG2_shWTAP_m6a.xls ../xls/HepG2_shCont_m6a.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$13=log(1/$13)/log(2);print $1,$2,$3,$13;}' OFS="\t" temp.xls > HepG2_m6a_shWTAP_FC_circos.txt

awk '{gsub(/chr/,"hs",$1);print $1,$2,$3,$15;}' OFS="\t" ../xls/HepG2_shCont_m6a.xls  > HepG2_shCont_m6a_circos.txt

rm -f temp.xls

partial_ideogram_conf="
<ideogram>\n
\n
<spacing>\n
default = 0.002r\n
break   = 0.2r\n
</spacing>\n
\n
<<include ideogram.position.conf>>\n
<<include ideogram.label.conf>>\n
<<include bands.conf>>\n
\n
radius*       = 0.92r\n
\n
</ideogram>\n

"
partial_ideogram_label_conf="
show_label       = yes\n
label_font       = bold\n
label_radius     = dims(image,radius)-95p\n
label_size       = 80\n
label_parallel   = yes\n
label_case       = lower\n
label_format     = eval(sprintf(\"chr%s\",var(label)))\n
"
all_ideogram_conf="
<ideogram>\n
\n
<spacing>\n
default = 0.002r\n
break   = 0.2r\n
</spacing>\n
\n
<<include ideogram.position.conf>>\n
<<include ideogram.label.conf>>\n
<<include bands.conf>>\n
\n
radius*       = 0.92r\n
\n
</ideogram>\n

"
all_ideogram_label_conf="
show_label       = yes\n
label_font       = bold\n
label_radius     = dims(image,radius)-95p\n
label_size       = 80\n
label_parallel   = yes\n
label_case       = upper\n
label_format     = eval(sprintf(\"%s\",var(label)))\n
"


cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/circos/conf/histogram/
echo -e $partial_ideogram_conf > ideogram.conf
echo -e $partial_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.partial.FC.conf > /dev/null 2>&1
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.FC.conf > /dev/null 2>&1
rm -f histogram.partial.FC.svg histogram.FC.svg

#cd /data/zhoukr/hhl_setd2_m6a/analysis/m6a-seq/HepG2/circos/conf/scatter/
#circos -conf scatter.partial.shCont_vs_shSetD2.conf > /dev/null 2>&1
#circos -conf scatter.shCont_vs_shSetD2.conf > /dev/null 2>&1
#rm -f scatter.partial.shCont_vs_shSetD2.svg scatter.shCont_vs_shSetD2.svg

