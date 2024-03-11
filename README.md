# 2021 Untreated JDM PBMCs, skin, and muscle

This project used PBMCs from controls, individuals with Juvenile Dermatomyositis (JDM), and some of the same patients when they were clinically inactive. A second cohort involved skin and muscle samples from children with JDM and controls without autoimmune disease. This repository contains the code to download the associated data on FigShare (including individual gene count data) and run the analysis in R Studio.

[Preprint](https://www.biorxiv.org/content/10.1101/2021.05.07.443007v2)

Paper: [Transcriptomes of peripheral blood mononuclear cells from juvenile dermatomyositis patients show elevated inflammation even when clinically inactive](https://pubmed.ncbi.nlm.nih.gov/34997119/)

## Data

[FigShare project](https://figshare.com/projects/2019_Untreated_juvenile_dermatomyositis_JDM_RNA-Seq/63539)

## Requirements - R libraries
* cowplot
* DESeq2
* ggalluvial
* ggbeeswarm
* ggrepel
* here
* openxlsx
* reshape2
* tidyverse
* UpSetR
* WGCNA

## Running analysis
1. Install R and R Studio if needed.
2. Clone this repository.
3. Install required R packages.
4. Knit the R Markdown files in order.
