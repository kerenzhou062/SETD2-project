#!/bin/sh

cd /data/zhoukr/hhl_setd2_m6a/mouse_ESCs_ChIP-seq/
## xls
rm -rf /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/xls
mkdir /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/xls
macsXlsPeakFilter.pl -pval 0.05 -original -x Ctrl_D0/macs1.4-merge/shCont_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/xls/mESCs_shCont_chip.xls
macsXlsPeakFilter.pl -pval 0.05 -original -x SetD2-KD_D0/macs1.4-merge/shSetD2_peaks.xls -o /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/xls/mESCs_shSetD2_chip.xls

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/circos/
macs2PeaksFC.pl -nonOverlapA -nonOverlapB -x ../xls/mESCs_shSetD2_chip.xls ../xls/mESCs_shCont_chip.xls -o temp.xls
awk '{gsub(/chr/,"mm",$1);$7=log(1/$7)/log(2);print $1,$2,$3,$7;}' OFS="\t" temp.xls > mESCs_chip_shSetD2_FC_circos.txt
rm -f temp.xls

awk '{gsub(/chr/,"mm",$1);print $1,$2,$3,$8;}' OFS="\t" ../xls/mESCs_shCont_chip.xls > mESCs_macs_shCont_circos.txt
awk '{gsub(/chr/,"mm",$1);print $1,$2,$3,$8;}' OFS="\t" ../xls/mESCs_shSetD2_chip.xls > mESCs_macs_shSetD2_circos.txt

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

cd /data/zhoukr/hhl_setd2_m6a/analysis/chip-seq/mESCs/circos/conf/heatmap/
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf heatmap.shCont_vs_shSetD2.conf > /dev/null 2>&1
rm -f heatmap.shCont_vs_shSetD2.svg
