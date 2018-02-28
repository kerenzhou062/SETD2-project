## [Download software](http://software.broadinstitute.org/gsea/downloads.jsp) ##
```
gsea2-2.2.4.jar
```

## Preparing gene_set ##
```
download from GSEA
```

## Preparing datasets ##
- [.cls](http://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats) files 
```
4	4	1
#shContD0	shContD6	shSetD2D0	shSetD2D6
shContD0	shContD6	shSetD2D0	shSetD2D6
```
- .txt expression files
```
NAME	DESCRIPTION	shContD0	shContD6	shSetD2D0	shSetD2D6
GNAI3	na	15.37	24.44	16.44	12.18
CDC45	na	17.74	13.27	18.81	11.62
...
```

## Running GSEA analysis ##

### Required fields ###
- Number of permutations `1000`
- Clappse dataset to gene symbols `false`
- Permutation type `phenotype`
- Chip platform(s) `no platform used`

### Basic fields ###
- Analysis name
```
gene_level_stem_cell_related_shSetD2D6_vs_shContD6
```
- Enrichment statistic `weighted`
- Metric for ranking genes `Ratio_of_Classes`
- Gene list sorting mode `real`
- Gene list ordering mode `descending`
- Max size:exclude larger sets `500`
- Min aiw:exclude smaller sets `1`


