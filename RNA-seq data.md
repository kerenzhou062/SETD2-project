The data processing procedure for RNA-seq data (Input samples of m6A-seq) were listed as below:

---
## Estimate exression level of genes  ##
We use RSEM (version 1.2.31) to estimate gene expression levels.
```bash
rsem-calculate-expression  --phred33-quals \
                           --strandedness reverse \
                           --bowtie2 \
                           --output-genome-bam \
                           -p 8 \
                           /data/zhoukr/hhl_setd2_m6a/HepG2_m6A-seq/HepG2_m6A-seq_shCont_input_rep1.fastq \
                           /data/zhoukr/hhl_setd2_m6a/reference/rsem-index/bowtie2-hg19/hg19 \
                           /data/zhoukr/hhl_setd2_m6a/analysis/RNA-seq/HepG2/RSEM/shCont_input_rep1

```

## Calculate differenetially expressive genes ##
We used featureCounts and DEGSeq to calculate DE genes.

```bash
featureCounts -T 10 -a gencode.v24lift37.annotation.gtf -g gene_id -F GTF -t gene -M -o temp.txt \
  input1.bam input2.bam ...
   > /dev/null 2>&1

#related bash scripts
*DEGseq.sh
```

