## Hickey et al. adolescent vaginal microbiome study

This repository contains the data and code necessary to reproduce the analyses and figures presented in the paper "Vaginal microbiota of adolescent girls resemble those of reproductive-age women prior to the onset of menarche" by Hickey et al. in preparation for submission (October 2014).

The analyses are separated into five sections available in regular Markdown (.md) and R Markdown (.Rmd) format:
* 01-data-prep includes the initial setup of taxon data, participant metadata, and custom functions and color palettes.
* 02-hclust-pcoa covers the hierarchical clustering and ordination (PCoA) approaches used to compare vaginal microbiota samples.
* 03-community-dynamics generates Appendix S2 and goes through exploratory data analysis of vaginal microbiota dynamics.
* 04-lmm-lab-ph includes the linear mixed-effects modeling of trends in lactic acid bacteria proportions and pH in relation to participant metadata.
* 05-vagina-vulva-comparison includes analyses performed to compare vaginal and vulvar microbiota samples.

View the Markdown files to see the full code and embedded figures.

If you want to run the analyses yourself, you can fork or clone this repository and run the RMarkdown scripts from the main directory. Data are in the data subdirectory, smaller scripts are in the scripts subdirectory, and output from my own run-through are stored in the various other subdirectories.

Please direct any inquiries regarding the analyses to Roxana Hickey at <roxana.hickey@gmail.com>.
