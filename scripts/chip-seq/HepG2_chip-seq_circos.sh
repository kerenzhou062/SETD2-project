#!/bin/sh
cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/circos/

macs2PeaksFC.pl -cutoff 0.1 -nonOverlapA -nonOverlapB -x ../xls/HepG2_shSetD2_chip.xls ../xls/HepG2_shCont_chip.xls -o temp.xls
awk '{gsub(/chr/,"hs",$1);$7=log(1/$7)/log(2);print $1,$2,$3,$7;}' OFS="\t" temp.xls > HepG2_chip_shSetD2_FC_circos.txt
rm -f temp.xls

awk '{gsub(/chr/,"hs",$1);print $1,$2,$3,$8;}' OFS="\t" ../xls/HepG2_shCont_chip.xls > HepG2_macs_shCont_circos.txt
awk '{gsub(/chr/,"hs",$1);print $1,$2,$3,$8;}' OFS="\t" ../xls/HepG2_shSetD2_chip.xls > HepG2_macs_shSetD2_circos.txt

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

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/HepG2/circos/conf/heatmap/
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf heatmap.shCont_vs_shSetD2.conf > /dev/null 2>&1
rm -f heatmap.shCont_vs_shSetD2.svg

