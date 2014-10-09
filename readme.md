## Hickey et al. adolescent vaginal microbiome study

This repository contains the data and code necessary to reproduce the analyses and figures presented in the paper "Vaginal microbiota of adolescent girls resemble those of reproductive-age women prior to the onset of menarche" by Hickey et al. in review as of October 2014.

The analyses are separated into five files available as both regular Markdown (.md) and R Markdown (.Rmd) files:
* 01-data-prep.Rmd includes the initial setup of taxon data, participant metadata, and custom functions and color palettes.
* 02-hclust-pcoa.Rmd covers the hierarchical clustering and ordination (PCoA) approaches used to compare vaginal microbiota samples.
* 03-community-dynamics.Rmd generates Appendix S2 and goes through exploratory data analysis of vaginal microbiota dynamics.
* 04-lmm-lab-ph.Rmd includes the linear mixed-effects modeling of trends in lactic acid bacteria pro- portions and pH in relation to participant metadata.
* 05-vagina-vulva-comparison.Rmd includes analyses performed to compare vaginal and vulvar mi- crobiota samples.

These are also published with embedded figures on RPubs at the following links, corresponding to the above Rmd files. These are useful if you are interested in following the analyses in detail without running the scripts yourself.
* http://rpubs.com/roxanahickey/adolescent-01-data-prep
* http://rpubs.com/roxanahickey/adolescent-02-hclust-pcoa
* http://rpubs.com/roxanahickey/adolescent-03-community-dynamics
* http://rpubs.com/roxanahickey/adolescent-04-lmm-lab-ph
* http://rpubs.com/roxanahickey/adolescent-05-vagina-vulva-comparison

If you do want to run the analyses yourself, everything you need is in this repository. Data are in the data subdirectory, smaller scripts are in the scripts subdirectory, and output from my own running of the analyses are stored in the various other subdirectories.

Please direct any inquiries regarding the analyses to Roxana Hickey at <roxana.hickey@gmail.com>.