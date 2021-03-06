## Vaginal and vulvar microbiota of perimenarcheal girls

![alt text](https://roxanahickey.files.wordpress.com/2014/10/silhouettes-lifetime-yellow2.png)

This repository contains the data and code necessary to reproduce the analyses and figures presented in the following publication:

Hickey RJ, Zhou X, Settles ML, Erb J, Malone K, Hansmann MA, Shew ML, Van Der Pol B, Fortenberry JD, Forney LJ. 2015. Vaginal microbiota of adolescent girls prior to the onset of menarche resemble those of reproductive-age women. mBio 6(2):e00097-15. [doi:10.1128/mBio.00097-15.](http://mbio.asm.org/content/6/2/e00097-15.long)

This work was presented previously in [talk](http://www.slideshare.net/roxana_hickey/hickeyuometa2014talk) and [poster](http://www.slideshare.net/roxana_hickey/hickey-isme15-poster) format.

Raw 454 sequence data are available via NCBI SRA/BioProject:
* BioProject accession [PRJNA266340](http://www.ncbi.nlm.nih.gov/bioproject/PRJNA266340)
* SRA accession [SRP055743](http://www.ncbi.nlm.nih.gov/Traces/sra/?study=SRP055743)

The analyses are separated into six sections available in regular Markdown (.md) and R Markdown (.Rmd) format (click on the links to view the Markdown files to see the full code and embedded figures):
* [01-data-prep](https://github.com/roxanahickey/adolescent/blob/master/01-data-prep.md) includes the initial setup of taxon data, participant metadata, and custom functions and color palettes.
* [02-hclust-pcoa](https://github.com/roxanahickey/adolescent/blob/master/02-hclust-pcoa.md) covers the hierarchical clustering and ordination (PCoA) approaches used to compare vaginal microbiota samples.
* [03-community-dynamics](https://github.com/roxanahickey/adolescent/blob/master/03-community-dynamics.md) generates community composition/metadata profiles and goes through exploratory data analysis of vaginal microbiota dynamics.
* [04-lmm-lab-ph](https://github.com/roxanahickey/adolescent/blob/master/04-lmm-lab-ph.md) includes the linear mixed-effects modeling of trends in lactic acid bacteria proportions and pH in relation to participant metadata.
* [05-vagina-vulva-comparison](https://github.com/roxanahickey/adolescent/blob/master/05-vagina-vulva-comparison.md) includes analyses performed to compare vaginal and vulvar microbiota samples.
* [06-post-review](https://github.com/roxanahickey/adolescent/blob/master/06-post-review.md) includes additional analyses performed following an initial round of peer review. (added 2015-01-12)

If you want to run the analyses yourself, you can fork or clone this repository and run the R Markdown scripts from the main directory. Please direct any inquiries regarding the analyses to Roxana Hickey at <roxana.hickey@gmail.com>.
