## Vaginal and vulvar microbiota of perimenarcheal girls

![alt text](https://roxanahickey.files.wordpress.com/2014/10/silhouettes-lifetime-yellow2.png)

This repository contains the data and code necessary to reproduce the analyses and figures presented in the paper "Vaginal microbiota of adolescent girls resemble those of reproductive-age women prior to the onset of menarche" by Hickey et al. in preparation for submission (October 2014).

This work has been presented previously in [talk](http://www.slideshare.net/roxana_hickey/hickeyuometa2014talk) and [poster](http://www.slideshare.net/roxana_hickey/hickey-isme15-poster) format.

The analyses are separated into five sections available in regular Markdown (.md) and R Markdown (.Rmd) format (click on the links to view the Markdown files to see the full code and embedded figures):
* [01-data-prep](https://github.com/roxanahickey/adolescent/blob/master/01-data-prep.md) includes the initial setup of taxon data, participant metadata, and custom functions and color palettes.
* [02-hclust-pcoa](https://github.com/roxanahickey/adolescent/blob/master/02-hclust-pcoa.md) covers the hierarchical clustering and ordination (PCoA) approaches used to compare vaginal microbiota samples.
* [03-community-dynamics](https://github.com/roxanahickey/adolescent/blob/master/03-community-dynamics.md) generates Appendix S2 and goes through exploratory data analysis of vaginal microbiota dynamics.
* [04-lmm-lab-ph](https://github.com/roxanahickey/adolescent/blob/master/04-lmm-lab-ph.md) includes the linear mixed-effects modeling of trends in lactic acid bacteria proportions and pH in relation to participant metadata.
* [05-vagina-vulva-comparison](https://github.com/roxanahickey/adolescent/blob/master/05-vagina-vulva-comparison.md) includes analyses performed to compare vaginal and vulvar microbiota samples.

If you want to run the analyses yourself, you can fork or clone this repository and run the R Markdown scripts from the main directory. Data are in the data subdirectory, shorter scripts are in the scripts subdirectory, and output from my own run-through of the analyses are stored in the various other subdirectories.

Please direct any inquiries regarding the analyses to Roxana Hickey at <roxana.hickey@gmail.com>.
