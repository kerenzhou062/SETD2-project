#!/bin/sh

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

cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/HepG2/H3K36me3-m6a/histogram
echo -e $partial_ideogram_conf > ideogram.conf
echo -e $partial_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.partial.FC.conf > /dev/null 2>&1
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.FC.conf > /dev/null 2>&1
rm -f histogram.partial.FC.svg histogram.FC.svg

#cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/HepG2/H3K36me3-m6a/scatter
#circos -conf scatter.partial.shCont_vs_shSetD2.conf > /dev/null 2>&1
#circos -conf scatter.shCont_vs_shSetD2.conf > /dev/null 2>&1
#rm -f scatter.partial.shCont_vs_shSetD2.svg scatter.shCont_vs_shSetD2.svg



cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/Hela/H3K36me3-m6a/histogram
echo -e $partial_ideogram_conf > ideogram.conf
echo -e $partial_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.partial.FC.conf > /dev/null 2>&1
echo -e $all_ideogram_conf > ideogram.conf
echo -e $all_ideogram_label_conf > ideogram.label.conf
circos -conf histogram.FC.conf > /dev/null 2>&1
rm -f histogram.partial.FC.svg histogram.FC.svg

#cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/Hela/H3K36me3-m6a/scatter
#circos -conf scatter.partial.shCont_vs_shSetD2.conf > /dev/null 2>&1
#circos -conf scatter.shCont_vs_shSetD2.conf > /dev/null 2>&1
#rm -f scatter.partial.shCont_vs_shSetD2.svg scatter.shCont_vs_shSetD2.svg

cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/HepG2/H3K36me3-m6a/heatmap
circos -conf heatmap.m6A_H3K36me3.conf > /dev/null 2>&1
rm -f heatmap.m6A_H3K36me3.svg

cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/HepG2/H3K36me3-m6a/heatmap
circos -conf heatmap.m6A_H3K36me3.conf > /dev/null 2>&1
rm -f heatmap.m6A_H3K36me3.svg

circos -conf heatmap.m6A_H3K36me3_partial.conf > /dev/null 2>&1
rm -f heatmap.m6A_H3K36me3_partial.svg

cd /data/zhoukr/hhl_setd2_m6a/analysis/circos/HepG2/H3K9me3-H3K36me3-m6a/heatmap
circos -conf heatmap.m6A_H3K36me3_H3K9me3.conf > /dev/null 2>&1
rm -f heatmap.m6A_H3K36me3_H3K9me3.svg

circos -conf heatmap.m6A_H3K36me3_H3K9me3_partial.conf > /dev/null 2>&1
rm -f heatmap.m6A_H3K36me3_H3K9me3_partial.svg
