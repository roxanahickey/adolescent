# Supplemental File 3. Analysis of vaginal microbiota dynamics in relation to pubertal development
Roxana J. Hickey  

# Description
This is a supplement to the paper "Vaginal microbiota of adolescent girls resemble those of reproductive-age women prior to the onset of menarche" by Hickey et al. Please refer to the paper for complete information about the objectives and study design.

The embedded R code works through the first set of analyses and generation of figures related to the assessment of vaginal microbiota composition in girls and mothers. The analyses can be run directly from the R Markdown file "adolescent-supp-03.Rmd" using RStudio. It should be run after "adolescent-supp-01.Rmd" and "adolescent-supp-02.Rmd".

## Objective
The first set of analyses dealt with characterizing vaginal microbiota composition and identifying major groups using clustering and ordination approaches. Next, we will take a closer look at the dynamics of these communities over time as girls progress through puberty and menarche. Special attention is devoted to assessing trends in the relative abundance of lactic acid bacteria and vaginal pH, as these are generally considered indicators of a 'healthy' vaginal microbiota in reproductive age women.

***
# Initial setup

Clear the workspace, load data from the previous step, and load necessary packages.

*Note: If you run the R Markdown script 'as is' from the same directory containing it and the 'data-input' and 'data-postproc' subdirectories, all figures will be printed inside the resulting PDF or HTML output. If you want to save the figures as individual files, uncomment any lines throughout the script starting with 'ggsave()' or 'pdf()'. I made note of each of these within the chunk code.*


```r
## Clear current workspace
rm(list=ls())

## Load the RData file created from adolescent-supp-02.Rmd
## (this contains data from supp-01 as well)
load("data-postproc/02-cluster-pcoa-last-run.RData")

## Load packages
library(ape)
library(ggplot2)
library(gplots)
```

```
## KernSmooth 2.23 loaded
## Copyright M. P. Wand 1997-2009
## 
## Attaching package: 'gplots'
## 
## The following object is masked from 'package:stats':
## 
##     lowess
```

```r
library(grid)
library(gridExtra)
library(reshape)
library(scales)

## Write new directories for output
dir.create("community-dynamics")
```

```
## Warning: 'community-dynamics' already exists
```

```r
dir.create("regression-analysis")
```

```
## Warning: 'regression-analysis' already exists
```

***
# Part I: Summarize community dynamics of vaginal (and vulvar) microbiota

First we will create an individual summary for each individual participant (and her mother, if applicable) to look at community composition of the vaginal and vulvar microbiota, along with selected metadata, at every visit during the study. The output is a PDF file with each page containing the following plots for a given subject:

* Community composition of girl's vaginal microbiota
* Community composition of girl's vulvar microbiota
* Community composition of mother's vaginal microbiota
* Girl's menarcheal status
* Girl's Tanner breast and genital stages (clinician-assessed)
* Vaginal pH of girl and mother (when available)

## Supplemental File 4. Individual summaries of vaginal and vulvar microbiota composition and associated metadata for all study participants


```r
## g_legend function from http://stackoverflow.com/questions/11883844/
## inserting-a-table-under-the-legend-in-a-ggplot2-histogram
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

pdf(file="community-dynamics/adolescent-supp-04.pdf", 
    width=8.5, height=11, pointsize=8)
for (i in as.numeric(unique(meta$gm.pair))){
  pick <- which(as.numeric(meta$gm.pair) == i)
  if (length(pick) >= 1){
    prop.tmp <- prop.red[,pick]
    meta.tmp <- meta[pick,]
    
    meta.tmp$race <- gsub("bla", "Black", meta.tmp$race)
    meta.tmp$race <- gsub("cauc", "Caucasian", meta.tmp$race)
    meta.tmp$race <- gsub("ami", "Native American", meta.tmp$race)
    meta.tmp$ethn <- gsub("nhis", "non-Hispanic", meta.tmp$ethn)
    meta.tmp$ethn <- gsub("hisp", "Hispanic", meta.tmp$ethn)
    
    taxa.pick25 <- order(rowMeans(as.matrix(prop.tmp), 
                                  na.rm=TRUE), 
                         decreasing=T)[1:25]
    prop.top25.tmp <- prop.tmp[taxa.pick25,]
    
    meta.prop25.tmp <- data.frame(cbind(meta.tmp$type.site, 
                                        meta.tmp$visit, 
                                        meta.tmp$men.stat,
                                        meta.tmp$tan.br.dr,
                                        meta.tmp$tan.gen.dr,
                                        meta.tmp$ph,
                                        t(prop.top25.tmp)))
    colnames(meta.prop25.tmp)[1:6] <- c("type.site", 
                                        "visit",
                                        "men.stat",
                                        "tan.br.dr",
                                        "tan.gen.dr",
                                        "ph")
    
    meta.prop25.lg.tmp <- melt(meta.prop25.tmp,
                               id.vars=c("type.site", 
                                         "visit",
                                         "men.stat",
                                         "tan.br.dr",
                                         "tan.gen.dr",
                                         "ph"))
    colnames(meta.prop25.lg.tmp)[7:8] <- c("taxon", "prop")
    
    meta.prop25.lg.tmp$visit <- as.integer(unfactor(meta.prop25.lg.tmp$visit))
    meta.prop25.lg.tmp$prop <- as.numeric(unfactor(meta.prop25.lg.tmp$prop))
    
    if (length(levels(meta.prop25.lg.tmp$type.site)) == 2){
      levels(meta.prop25.lg.tmp$type.site) <- c("Girl Vagina", 
                                                "Girl Vulva")}
    else if (length(levels(meta.prop25.lg.tmp$type.site)) == 3){
      levels(meta.prop25.lg.tmp$type.site) <- c("Girl Vagina", 
                                                "Girl Vulva", 
                                                "Mother Vagina")}
    
    ## Community composition
    gg.taxa <- ggplot(meta.prop25.lg.tmp, 
                      aes(x=visit, y=prop, fill=taxon))
    gg.taxa.bar <- gg.taxa + geom_bar(stat="identity") +
      facet_wrap(~ type.site, ncol=1, scale="fixed") +
      scale_fill_manual(values=col.taxa[taxa.pick25], name="Taxon") +
      scale_x_discrete(limits=seq(1, max(meta.prop25.lg.tmp$visit), by=1), 
                       name="") +
      ylab("Proportion") +
      ggtitle(paste(unique(meta.tmp$gm.pair.fullID), ": ",
                    unique(meta.tmp$race), ", ",
                    unique(meta.tmp$ethn), ", ", 
                    round(min(na.exclude(meta.tmp$age.sampling)), digits=1),
                    " y/o at enrollment\n", sep="")) +
      theme_cust +
      theme(axis.text.x=element_blank(),
            legend.text=element_text(size=6),
            legend.title=element_text(size=8),
            legend.position=c(0,0.5),
            legend.justification=c(0,0.5))
    
    ## Now we need to go back to the original metadata table which has each
    ## visit listed in its own row. This way we don't have duplicate menarche
    ## status and Tanner scores to plot (one matching vagina, one matching
    ## vulva when they are identical)
    meta2.tmp <- subset(meta.orig, 
                        sample.ID.vag %in% rownames(meta.tmp) | 
                          sample.ID.vul %in% rownames(meta.tmp))
    
    meta2.tmp$men.stat <- factor(meta2.tmp$men.stat, 
                                  levels=c(NULL,"post","pre"))
    
    ## Menarche status
    gg.men.stat <- ggplot(subset(meta2.tmp, type=="girl"), 
                          aes(x=visit, y=subject, color=men.stat))
    gg.men.stat.dot <- gg.men.stat + 
      geom_point(size=5) +  
      scale_color_manual(values=col.men.stat[2:3],
                         breaks=c("pre","post"),
                         labels=c("Pre","Post"),
                         na.value="gray70",
                         name="Menarche status",
                         drop=FALSE) +
      scale_x_discrete(limits=seq(1, max(meta.prop25.lg.tmp$visit), by=1), 
                       name="") +
      scale_y_discrete(labels="Menarche", name="") +
      theme_cust +
      theme(axis.text.x=element_blank(),
            legend.text=element_text(size=6),
            legend.title=element_text(size=8),
            legend.position=c(0,1),
            legend.justification=c(0,1))
    
    ## Tanner
    tanner.tmp <- meta2.tmp[,c("visit",
                               "tan.br.dr",
                               "tan.gen.dr",
                               "age.sampling")]
    
    tanner.lg.tmp <- melt(tanner.tmp, id.vars=c("visit", "age.sampling"))
    
    tanner.lg.tmp$value <- factor(tanner.lg.tmp$value, 
                                  levels=c("1","2","3","4","5"))
    
    gg.tanner <- ggplot(tanner.lg.tmp, 
                        aes(x=visit, y=variable, color=value))
    gg.tanner.dot <- gg.tanner +  
      geom_point(size=5) +  
      scale_color_manual(values=col.tanner,
                         name="Tanner stage",
                         drop=FALSE) +  
      scale_x_discrete(limits=seq(1, max(tanner.lg.tmp$visit), by=1),
                       name="Visit") +
      scale_y_discrete(breaks=c("tan.gen.dr","tan.br.dr"), 
                       labels=c("Genital", "Breast"), 
                       name="") +
      theme_cust +
      theme(legend.text=element_text(size=6),
            legend.title=element_text(size=8),
            legend.position=c(0,1),
            legend.justification=c(0,1))
    
    # pH
    gg.ph <- ggplot(subset(meta2.tmp, type=="girl"), 
                    aes(x=visit, y=ph))
    gg.ph.line <- gg.ph +
      geom_point(size=2) +
      geom_line(lty=3) +
      scale_x_discrete(limits=seq(1, max(meta2.tmp$visit), by=1), 
                       name="Visit") +
      scale_y_continuous(limits=c(4.0, 8.0)) +
      geom_text(label=meta2.tmp$ph[meta2.tmp$type=="girl"], 
                vjust=-0.5, size=3.5) +
      ylab("Vaginal pH") +
      theme_cust_nominor
    
    check.ph <- is.na(meta2.tmp$ph[meta2.tmp$type=="girl"]) 

    if  ("FALSE" %in% check.ph) {
      grid.arrange(gg.taxa.bar + theme(legend.position="none",
                                       plot.margin=unit(c(0.25, 0.5, 0.25, 2), "lines")),
                   g_legend(gg.taxa.bar), 
                   gg.men.stat.dot + theme(legend.position="none",
                                           plot.margin=unit(c(0.25, 0.5, 0.25, 0), "lines")),
                   g_legend(gg.men.stat.dot),
                   gg.tanner.dot + theme(legend.position="none",
                                         axis.text.x=element_blank(),
                                         axis.title.x=element_blank(),
                                         plot.margin=unit(c(0.25, 0.5, 0.25, 1), "lines")),
                   g_legend(gg.tanner.dot),
                   gg.ph.line + theme(legend.position="none",
                                      plot.margin=unit(c(0.25, 0.5, 0.25, 3.25), "lines")),
                   ncol=2, widths=c(4,1), heights=c(6, 1, 1.25, 1.25))
      }
    else 
      grid.arrange(gg.taxa.bar + theme(legend.position="none",
                                       plot.margin=unit(c(0.25, 0.5, 0.25, 2), "lines")),
                   g_legend(gg.taxa.bar), 
                   gg.men.stat.dot + theme(legend.position="none",
                                           plot.margin=unit(c(0.25, 0.5, 0.25, 0), "lines")),
                   g_legend(gg.men.stat.dot),
                   gg.tanner.dot + theme(legend.position="none",
                                         plot.margin=unit(c(0.25, 0.5, 0.25, 1), "lines")),
                   g_legend(gg.tanner.dot),
                   ncol=2, widths=c(4,1), heights=c(7, 1, 1.5))
    }
  }
```

```
## Warning: Removed 10 rows containing missing values (geom_point).
## Warning: Removed 7 rows containing missing values (geom_path).
## Warning: Removed 10 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 4 rows containing missing values (geom_point).
## Warning: Removed 3 rows containing missing values (geom_path).
## Warning: Removed 4 rows containing missing values (geom_text).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_text).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_text).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_text).
## Warning: Removed 5 rows containing missing values (geom_point).
## Warning: Removed 5 rows containing missing values (geom_path).
## Warning: Removed 5 rows containing missing values (geom_text).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 2 rows containing missing values (geom_text).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_path).
## Warning: Removed 2 rows containing missing values (geom_text).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_path).
## Warning: Removed 2 rows containing missing values (geom_text).
## Warning: Removed 6 rows containing missing values (geom_point).
## Warning: Removed 5 rows containing missing values (geom_path).
## Warning: Removed 6 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 6 rows containing missing values (geom_point).
## Warning: Removed 6 rows containing missing values (geom_path).
## Warning: Removed 6 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_path).
## Warning: Removed 1 rows containing missing values (geom_text).
## Warning: Removed 2 rows containing missing values (geom_point).
## Warning: Removed 2 rows containing missing values (geom_path).
## Warning: Removed 2 rows containing missing values (geom_text).
## Warning: Removed 4 rows containing missing values (geom_point).
## Warning: Removed 4 rows containing missing values (geom_path).
## Warning: Removed 4 rows containing missing values (geom_text).
## Warning: Removed 3 rows containing missing values (geom_point).
## Warning: Removed 3 rows containing missing values (geom_path).
## Warning: Removed 3 rows containing missing values (geom_text).
## Warning: Removed 1 rows containing missing values (geom_point).
## Warning: Removed 1 rows containing missing values (geom_path).
## Warning: Removed 1 rows containing missing values (geom_text).
```

```r
dev.off()
```

```
## pdf 
##   2
```

## Highlight examples of _Lactobacillus_ dominant vaginal microbiota

Next we will plot just a few examples of vaginal microbiota with _Lactobacillus_ predominant by menarche. I have selected subjects 102, 103, 107 and 109. We will plot only the vaginal microbiota composition, menarche status, and Tanner breast stage for each. This also requires a bit of setup:


```r
## Subset prop and meta for subjects 102, 103, 107 and 109
meta.fig5 <- subset(meta, subject %in% c(102, 103, 107, 109) & site=="vag")
prop.fig5 <- prop.red[,colnames(prop.red) %in% rownames(meta.fig5)]

taxa.pick40 <- order(rowMeans(as.matrix(prop.fig5), na.rm=TRUE),
                     decreasing=T)[1:40]
prop.fig5.top40 <- prop.fig5[taxa.pick40,]
    
mp.fig5 <- data.frame(cbind(as.character(meta.fig5$subject),
                            meta.fig5$visit, 
                            as.character(meta.fig5$men.stat),
                            meta.fig5$tan.br.dr,
                            t(prop.fig5.top40)))
colnames(mp.fig5)[1:4] <- c("subject",
                            "visit",
                            "men.stat",
                            "tan.br.dr")

mp.fig5.lg <- melt(mp.fig5, id.vars=c("subject",
                                      "visit",
                                      "men.stat",
                                      "tan.br.dr"))

mp.fig5.lg$visit <- as.integer(unfactor(mp.fig5.lg$visit))
mp.fig5.lg$value <- as.numeric(unfactor(mp.fig5.lg$value))

mp.fig5.lg$tan.br.dr <- factor(mp.fig5.lg$tan.br.dr, levels=c(1,2,3,4,5))

mp.s102 <- subset(mp.fig5.lg, subject==102)
mp.s103 <- subset(mp.fig5.lg, subject==103)
mp.s107 <- subset(mp.fig5.lg, subject==107)
mp.s109 <- subset(mp.fig5.lg, subject==109)

## Plot vaginal microbiota composition of subject 102 (w/ legend)
gg.taxa.bar.102 <- ggplot(mp.s102, aes(x=visit, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=col.taxa[taxa.pick40], name="Taxon") +
  scale_x_discrete(limits=seq(1, max(mp.s102$visit), by=1), 
                   name="") +
  ylab("Proportion") +
  ggtitle("Subject 102") +
  theme_cust +
  guides(fill=guide_legend(ncol=8)) +
  theme(axis.text.x=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.position=c(0,0),
        legend.justification=c(0,0),
        legend.direction="horizontal")

## Plot vaginal microbiota composition of subject 103
gg.taxa.bar.103 <- ggplot(mp.s103, aes(x=visit, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=col.taxa[taxa.pick40], name="Taxon") +
  scale_x_discrete(limits=seq(1, max(mp.s103$visit), by=1), 
                   name="") +
  ylab("Proportion") +
  ggtitle("Subject 103") +
  theme_cust +
  theme(axis.text.x=element_blank(),
        legend.position="none")

## Plot vaginal microbiota composition of subject 107
gg.taxa.bar.107 <- ggplot(mp.s107, aes(x=visit, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=col.taxa[taxa.pick40], name="Taxon") +
  scale_x_discrete(limits=seq(1, max(mp.s107$visit), by=1), 
                   name="") +
  ylab("Proportion") +
  ggtitle("Subject 107") +
  theme_cust +
  theme(axis.text.x=element_blank(),
       legend.position="none")

## Plot vaginal microbiota composition of subject 109
gg.taxa.bar.109 <- ggplot(mp.s109, aes(x=visit, y=value, fill=variable)) +
  geom_bar(stat="identity") +
  scale_fill_manual(values=col.taxa[taxa.pick40], name="Taxon") +
  scale_x_discrete(limits=seq(1, max(mp.s109$visit), by=1), 
                   name="") +
  ylab("Proportion") +
  ggtitle("Subject 109") +
  theme_cust +
  theme(axis.text.x=element_blank(),
        legend.position="none")

## Plot Tanner/menarche of subject 102 (w/ legend)
gg.tb.dot.102 <- ggplot(mp.s102, aes(x=visit, y=subject, 
                                     shape=men.stat, color=tan.br.dr)) + 
  geom_point(size=5) +  
  scale_color_manual(values=col.tanner,
                     name="Tanner breast",
                     drop=FALSE) +
  scale_shape_manual(breaks=c("pre", "post"),
                     labels=c("Pre", "Post"),
                     values=c(16,1),
                     name="Menarche status") +
  scale_x_discrete(limits=seq(1, max(mp.s102$visit), by=1), 
                   name="Visit") +
  ylab("") +
  theme_cust +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.position=c(0.7,0.5),
        legend.justification=c(0.5,0.5),
        legend.direction="horizontal")

## Plot Tanner/menarche of subject 103
gg.tb.dot.103 <- ggplot(mp.s103, aes(x=visit, y=subject, 
                                     shape=men.stat, color=tan.br.dr)) + 
  geom_point(size=5) +  
  scale_color_manual(values=col.tanner,
                     name="Tanner breast",
                     drop=FALSE) +
  scale_shape_manual(breaks=c("pre", "post"),
                     labels=c("Pre", "Post"),
                     values=c(16,1),
                     name="Menarche status") +
  scale_x_discrete(limits=seq(1, max(mp.s103$visit), by=1), 
                   name="Visit") +
  ylab("") +
  theme_cust +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position="none")

## Plot Tanner/menarche of subject 107
gg.tb.dot.107 <- ggplot(mp.s107, aes(x=visit, y=subject, 
                                     shape=men.stat, color=tan.br.dr)) + 
  geom_point(size=5) +  
  scale_color_manual(values=col.tanner,
                     name="Tanner breast",
                     drop=FALSE) +
  scale_shape_manual(breaks=c("pre", "post"),
                     labels=c("Pre", "Post"),
                     values=c(16,1),
                     name="Menarche status") +
  scale_x_discrete(limits=seq(1, max(mp.s107$visit), by=1), 
                   name="Visit") +
  ylab("") +
  theme_cust +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position="none")

## Plot Tanner/menarche of subject 109
gg.tb.dot.109 <- ggplot(mp.s109, aes(x=visit, y=subject, 
                                     shape=men.stat, color=tan.br.dr)) + 
  geom_point(size=5) +  
  scale_color_manual(values=col.tanner,
                     name="Tanner breast",
                     drop=FALSE) +
  scale_shape_manual(breaks=c("pre", "post"),
                     labels=c("Pre", "Post"),
                     values=c(16,1),
                     name="Menarche status") +
  scale_x_discrete(limits=seq(1, max(mp.s109$visit), by=1), 
                   name="Visit") +
  ylab("") +
  theme_cust +
  theme(axis.ticks.y=element_blank(),
        axis.text.y=element_blank(),
        legend.position="none")
```

## Figure 5. Transitions to _Lactobacillus_-dominant vaginal microbiota in four perimenarcheal girls.
Each panel shows the vaginal bacterial community profiles and associated pubertal development of four participants sampled longitudinally. The top bar plot of each panel shows the proportions of bacterial taxa present at each sampling, with colors of taxa indicated in the legend at bottom left. Below each bar plot the menarcheal status (M), Tanner breast (TB) and Tanner genital (TG) scores as assessed by a clinician are indicated following the color scheme shown in the legend at bottom right. The participant’s age in years at the time of sampling is indicated below the bars in each panel. Empty spaces in the plots indicate a skipped scheduled quarterly visit or a visit in which a vaginal swab or metadata were not collected.


```r
## Uncomment next line to save as a PDF
# pdf("community-dynamics/fig-5-community-dynamics.pdf", width=12, height=8, pointsize=8)
grid.arrange(gg.taxa.bar.107 + theme(plot.margin=unit(c(0.25, 1, 0, 0.8), "lines")),
             gg.taxa.bar.109 + theme(plot.margin=unit(c(0.25, 1, 0, 0.8), "lines")),
             gg.tb.dot.107 + theme(plot.margin=unit(c(0.25, 1, 0.5, 2.4), "lines")),
             gg.tb.dot.109 + theme(plot.margin=unit(c(0.25, 1, 0.5, 2.4), "lines")),
             gg.taxa.bar.102 + theme(legend.position="none",
                                     plot.margin=unit(c(0.25, 1, 0, 0.8), "lines")),
             gg.taxa.bar.103 + theme(plot.margin=unit(c(0.25, 1, 0, 0.8), "lines")),
             gg.tb.dot.102 + theme(legend.position="none",
                                   plot.margin=unit(c(0.25, 1, 0.5, 2.4), "lines")),
             gg.tb.dot.103 + theme(plot.margin=unit(c(0.25, 1, 0.5, 2.4), "lines")),
             g_legend(gg.taxa.bar.102),
             g_legend(gg.tb.dot.102),
             ncol=2, heights=c(3,0.75,3,0.75,1.25))
```

![plot of chunk fig-5-community-dynamics](./adolescent-supp-03_files/figure-html/fig-5-community-dynamics.png) 

```r
dev.off()
```

```
## null device 
##           1
```

## _Gardnerella_ in perimenarcheal vaginal microbiota

_Gardnerella vaginalis_ was surprisingly common among girls in our study. This is notable because the species (or genus) is commonly associated with bacterial vaginosis, and some have argued it is acquired through sexual contact. However, the girls in this study had no history of sexual activity and reported themselves in good health. Below we plot the changes in proportion of _Gardnerella_ over time in girls who had at least 5% _Gardnerella_ in her vaginal microbiota at any point in time during the study.


```r
meta.vag$Gvag <- spe.prop.vag[,"Gardnerella_vaginalis"]
meta.vag$Gardnerella <- rowSums(spe.prop.vag[,grep("Gardnerella", colnames(spe.prop.vag))])

meta.vag.gard.05 <- subset(meta.vag, type=="girl" & Gardnerella >= 0.05)
unique(meta.vag.gard.05$subject) # n=11
```

```
##  [1] 102 105 112 115 116 120 121 124 126 129 135
## 55 Levels: 101 102 103 104 105 106 107 108 109 110 111 112 113 114 ... 235
```

```r
## We want to plot Gardnerella over time in all of the girls who carried it at some point, 
## so we need to grab all the obervations for the girls listed above.
meta.vag.gard.05.fullobs <- meta.vag[meta.vag$subject %in% meta.vag.gard.05$subject,]
```

## Figure S3. Proportion of _Gardnerella_ over time in the vaginal microbiota of girls.
_Gardnerella_ was present in the vaginal microbiota at a proportion of 0.05 or greater at least once in 11/31 adolescent participants. Each subplot shows the proportion of _Gardnerella_ (encompassing sequence reads assigned to either the species level as _G. vaginalis_ or genus level as _Gardnerella_) in the vaginal microbiota of a single participant at each clinical visit. Open circles represent premenarcheal samples, and filled circles represent postmenarcheal samples. The x-axis indicates the clinical visit at which each sample was collected; visits occurred approximately every three months.


```r
gg.girl.gard <- ggplot(meta.vag.gard.05.fullobs, 
                       aes(x=visit, y=Gardnerella, group=subject, 
                           color=type, shape=men.stat))
gg.girl.gard + geom_line(lty=3, alpha=0.5) +
  facet_wrap( ~ gm.pair.fullID, ncol=4) +
  geom_point(size=3) +
  scale_color_manual(values="#01665E", guide=FALSE) +
  scale_shape_manual(values=c(16, 1), 
                     breaks=c("pre", "post"), 
                     labels=c("Pre", "Post"),
                     name="Menarche\nStatus", na.value=10) +
  xlab("Visit") +
  ylab("Proportion of Gardnerella") +
  scale_x_discrete(1:14, name="Visit No.") +
  ylim(0,1) +
  theme_cust_nominor +
  theme(axis.text=element_text(size=6),
        legend.justification=c(1,1), 
        legend.position=c(0.95,0.3))
```

```
## geom_path: Each group consist of only one observation. Do you need to adjust the group aesthetic?
```

![plot of chunk fig-s3-gardnerella-vag](./adolescent-supp-03_files/figure-html/fig-s3-gardnerella-vag.png) 

```r
## Uncomment next line to save as a PDF
# ggsave("community-dynamics/fig-s3-gardnerella-vag.pdf", width=8, height=4, units="in")
```

## _Gardnerella_ in vaginal and vulvar microbiota

I was also curious to see whether the same patterns were observed in the vulva. Not surprisingly, they are fairly similar, although subject 135 is an interesting exception with a might higher proportion of _Gardnerella_ on the vulva compared to the vagina.


```r
meta$Gvag <- t(prop.red["Gardnerella_vaginalis",])
meta$Gardnerella <- colSums(prop.red[grep("Gardnerella", rownames(prop.red)),])

meta.girl.gard.05 <- subset(meta, type=="girl" & Gardnerella >= 0.05)
unique(meta.girl.gard.05$subject) # n=11
```

```
##  [1] 102 105 112 115 116 120 121 124 126 129 135
## 55 Levels: 101 102 103 104 105 106 107 108 109 110 111 112 113 114 ... 235
```

```r
## Subjects who had Gardnerella in vagina
unique(meta.girl.gard.05$subject[meta.girl.gard.05$site=="vag"]) # n=11
```

```
##  [1] 102 105 112 115 116 120 121 124 126 129 135
## 55 Levels: 101 102 103 104 105 106 107 108 109 110 111 112 113 114 ... 235
```

```r
## Subjects who had Gardnerella on vulva
unique(meta.girl.gard.05$subject[meta.girl.gard.05$site=="vul"]) # n=10
```

```
##  [1] 102 105 112 115 116 121 124 126 129 135
## 55 Levels: 101 102 103 104 105 106 107 108 109 110 111 112 113 114 ... 235
```

```r
meta.girl.gard.05.fullobs <- meta[meta$subject %in% meta.girl.gard.05$subject,]
```

Plot vaginal and vulvar samples of girls with at least 5% _Gardnerella_:

```r
gg.girl.gard <- ggplot(meta.girl.gard.05.fullobs, 
                       aes(x=visit, y=Gardnerella, group=site, color=site, shape=men.stat))
gg.girl.gard + geom_line(lty=3, alpha=0.5) +
  facet_wrap( ~ gm.pair.fullID, ncol=4) +
  geom_point(size=3) +
  scale_color_manual(values=col.site, name="Site", 
                     breaks=c("vag","vul"), 
                     labels=c("Vagina", "Vulva")) +
  scale_shape_manual(values=c(16, 1), name="Menarche\nStatus", 
                     breaks=c("pre", "post"),
                     labels=c("Pre", "Post")) +
  xlab("Visit") +
  ylab("Proportion of Gardnerella") +
  scale_x_discrete(1:14, name="Visit No.") +
  ylim(0,1) +
  theme_cust_nominor +
  theme(axis.text=element_text(size=6))
```

```
## geom_path: Each group consist of only one observation. Do you need to adjust the group aesthetic?
```

![plot of chunk gardnerella-vag-vul](./adolescent-supp-03_files/figure-html/gardnerella-vag-vul.png) 

***
# Part II: Multiple linear regression analysis of lactic acid bacteria and vaginal pH

Now we want to extract some general trends from the data. We are going to focus on trends in lactic acid bacteria and vaginal pH. We will use a fairly simple multiple linear regression approach to identify variables significantly related to these changes. First we should simplify our variables and select only one of the Tanner measures to use as a proxy for pubertal development (we currently have four: breast self-assessement, breast clinician-assessment, genital self-assessment, genital clinician assessment). To do so we simply calculate the correlations between Tanner scores by body site (breast/genital) and self/clinician assessment to see how well they agree.

## Correlation of Tanner scores


```r
## Spearman
cor.sp <- cor(meta[meta$type=="girl" & meta$site=="vag", c("tan.gen.dr","tan.gen.self","tan.br.dr","tan.br.self")], 
              use="pairwise.complete.obs", method="spearman")
cor.sp
```

```
##              tan.gen.dr tan.gen.self tan.br.dr tan.br.self
## tan.gen.dr       1.0000       0.6834    0.7273      0.5678
## tan.gen.self     0.6834       1.0000    0.5942      0.7072
## tan.br.dr        0.7273       0.5942    1.0000      0.6336
## tan.br.self      0.5678       0.7072    0.6336      1.0000
```

## Figure S4. Correlations of Tanner breast and genital scores according to self and clinician assessment.
Spearman’s correlation coefficients were calculated to determine how well breast/genital scores agreed for both clinician- (upper left) and self-assessment (upper right).  Clinician- and self-assessments for Tanner breast (lower left) and Tanner genital (lower right) scores were also correlated. Correlation coefficients are reported in the header of each plot. Clinician-assessed Tanner breast scores were primarily used for analyses described in the paper.


```r
## Tanner genital vs. Tanner breast (clinician)
gg.tan.dr <- ggplot(meta[meta$type=="girl" & meta$site=="vag",], aes(x=tan.br.dr, y=tan.gen.dr)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.2), color="#8B1C62", size=3, alpha=0.7) +
  ggtitle("Clinician Assessment\n(Spearman's rho = 0.73)") +
  xlab("Tanner Breast") +
  ylab("Tanner Genital") +
  theme_cust_nominor

## Tanner genital vs. Tanner breast (self)
gg.tan.self <- ggplot(meta[meta$type=="girl" & meta$site=="vag",], aes(x=tan.br.self, y=tan.gen.self)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.2), color="#8B1C62", size=3, alpha=0.7) +
  ggtitle("Self Assessment\n(Spearman's rho = 0.71)") +
  xlab("Tanner Breast") +
  ylab("Tanner Genital") +
  theme_cust_nominor

## Tanner breast (self) vs. Tanner breast (dr)
gg.tb.dr.self <- ggplot(meta[meta$type=="girl" & meta$site=="vag",], aes(x=tan.br.dr, y=tan.br.self)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.2), color="#8B1C62", size=3, alpha=0.7) +
  ggtitle("Tanner Breast, Self vs. Clinician\n(Spearman's rho = 0.63)") +
  xlab("Clinician") +
  ylab("Self") +
  theme_cust_nominor

## Tanner genital (self) vs. Tanner genital (dr)
gg.tg.dr.self <- ggplot(meta[meta$type=="girl" & meta$site=="vag",], aes(x=tan.gen.dr, y=tan.gen.self)) + 
  geom_jitter(position=position_jitter(width=0.2, height=0.2), color="#8B1C62", size=3, alpha=0.7) +
  ggtitle("Tanner Genital, Self vs. Clinician\n(Spearman's rho = 0.68)") +
  xlab("Clinician") +
  ylab("Self") +
  theme_cust_nominor

## Uncomment the next line to save as a PDF
# pdf("fig-s4-tanner-correlations.pdf", width=10, height=10)
multiplot(gg.tan.dr, gg.tan.self, gg.tb.dr.self, gg.tg.dr.self, layout=matrix(c(1,2,3,4), ncol=2, byrow=TRUE))
```

```
## Warning: Removed 12 rows containing missing values (geom_point).
## Warning: Removed 9 rows containing missing values (geom_point).
## Warning: Removed 8 rows containing missing values (geom_point).
## Warning: Removed 12 rows containing missing values (geom_point).
```

![plot of chunk fig-s4-tanner-correlations](./adolescent-supp-03_files/figure-html/fig-s4-tanner-correlations.png) 

```r
dev.off()
```

```
## null device 
##           1
```

As the plots above show, the Tanner scores (breast and genital, self and clinician) are fairly well correlated. This is not too surprising since the Tanner scores for breast and genital were designed to capture simultaneous stages of development. One of the inclusion criteria for participation in the study was Tanner breast of at least stage 2, so it makes intuitive sense to rely the Tanner breast score as a proxy for pubertal development in our study. Additionally, the clinician-assessed Tanner breast score has more complete observations (192 out of 198 vagina samples) than the clinician assessed Tanner genital score (186/198). Therefore, we'll just use the Tanner breast score for the subsequent analyses.

## Multiple linear regression
Our goal here is to determine whether the proportion of lactic acid bacteria and vaginal pH are significantly associated with Tanner stage, age and/or menarche status. We'll set up several regression models (for each of the two response variables) and select a model **HOW**

*Note: In order to perform linear regression on proportional data (e.g., LAB proportions), the data should be approximately normally distributed. However, our proportions of LAB are heavily skewed with many close to 1 and quite a few close to zero, with far less in the middle. Here we subject the data to a logit transform, log(y/(1-y)), described by Warton and Hui as a preferred method for ecological non-binomial data such as proportions. However, because we have 0's and 1's in our data matrix, these would result in some -Inf and Inf values that cause problems for the linear model. An ad hoc solution to this is to add a nominal value epsilon to the numerator and denominator during the logit transform. Because our data are more skewed toward proportions close to 1 (e.g., check this with hist(meta.vagLAB)), we take as epsilon the difference between 1 and the largest proportion less than 1. This is defined below as eps. For discussion and justification of this approach, see Warton and Hui (2011) The arcsine is asinine: the analysis of proportions in ecology. Ecology 92:3–10.*

### Setup for modeling LAB

```r
## First we add the proportions of individual Lactobacillus spp., Lactobacillus genus, and lactic acid bacteria as additional variables to the metadata table
meta$Lactobacillus_crispatus <- t(prop.red["Lactobacillus_crispatus",])
meta$Lactobacillus_gasseri <- t(prop.red["Lactobacillus_gasseri",])
meta$Lactobacillus_iners <- t(prop.red["Lactobacillus_iners",])
meta$Lactobacillus_jensenii <- t(prop.red["Lactobacillus_jensenii",])
meta$Lactobacillus <- colSums(prop.red[grep("Lactobacillus", rownames(prop.red)),])
meta$LAB <- colSums(prop.red[c(6, 28, 43:49, 72:74),]) # 6=Aerococcus, 28=Facklamia, 43:49=Lactobacillus, 72:74=Streptococcus

## Define modified logit transform so we can specify epsilon adjustment
mod.logit <- function(x, eps) { log((x+eps)/(1-x+eps)) }

## Calculate epsilon:
## Subtract values of LAB from 1 
dif.lab <- 1-meta$LAB

## Take the smallest absolute difference greater than zero (~0.000075)
eps.lab <- min(dif.lab[dif.lab > 0])
eps.lab
```

```
## [1] 7.471e-05
```

```r
## Add logit-transformed LAB proportions to metadata
meta$LAB.logit <- mod.logit(meta$LAB, eps.lab)

## Subset metadata to only girl vagina samples and exclude Tanner breast = 1 (only 3 samples)
meta.girl.vag <- subset(meta, type=="girl" & site=="vag" & tan.br.dr %in% c(NA,2,3,4,5))

## Set these variables as factors so the linear models treat them as such
meta.girl.vag$tan.br.dr <- factor(meta.girl.vag$tan.br.dr, ordered=TRUE) # set as ordered factor
meta.girl.vag$men.stat <- factor(meta.girl.vag$men.stat, levels=c("pre","post"))
meta.girl.vag$subject <- factor(meta.girl.vag$subject)
```

## Model comparisons for LAB
Now that we're all set up, we start testing regression models. The first is a very simple model containing four variables (Tanner stage, age, menarche status and subject):

```r
## First, a very simple model with four variables and no interactions:
fit.lab.1 <- lm(LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject, meta.girl.vag)
anova(fit.lab.1) 
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##               Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr      3    528   176.0   24.68 4.0e-13 ***
## age.sampling   1     91    90.7   12.71 0.00048 ***
## men.stat       1     26    25.9    3.63 0.05845 .  
## subject       29   1132    39.0    5.47 9.8e-13 ***
## Residuals    156   1113     7.1                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(fit.lab.1)
```

```
## 
## Call:
## lm(formula = LAB.logit ~ tan.br.dr + age.sampling + men.stat + 
##     subject, data = meta.girl.vag)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
##  -8.43  -1.22   0.00   1.36   5.25 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)   -12.374      6.033   -2.05  0.04192 *  
## tan.br.dr.L     1.612      0.849    1.90  0.05947 .  
## tan.br.dr.Q    -1.148      0.492   -2.33  0.02088 *  
## tan.br.dr.C     0.919      0.401    2.29  0.02337 *  
## age.sampling    1.351      0.487    2.77  0.00628 ** 
## men.statpost   -1.658      0.692   -2.40  0.01780 *  
## subject102     -3.451      1.194   -2.89  0.00441 ** 
## subject103     -0.613      1.108   -0.55  0.58047    
## subject104      0.351      1.486    0.24  0.81338    
## subject105     -3.608      1.577   -2.29  0.02352 *  
## subject106      1.687      1.197    1.41  0.16079    
## subject107     -0.474      1.196   -0.40  0.69248    
## subject108      2.994      1.200    2.49  0.01368 *  
## subject109      2.218      1.220    1.82  0.07104 .  
## subject110      2.498      3.204    0.78  0.43671    
## subject111      2.724      1.249    2.18  0.03065 *  
## subject112      0.619      1.373    0.45  0.65267    
## subject113      1.265      1.422    0.89  0.37512    
## subject114      3.681      1.472    2.50  0.01343 *  
## subject115     -2.626      1.775   -1.48  0.14108    
## subject116     -1.867      1.775   -1.05  0.29438    
## subject118    -13.247      2.818   -4.70  5.7e-06 ***
## subject120     -4.910      1.444   -3.40  0.00085 ***
## subject121     -0.533      1.459   -0.37  0.71526    
## subject123      0.504      1.457    0.35  0.72993    
## subject124      3.280      1.623    2.02  0.04502 *  
## subject125      2.599      1.867    1.39  0.16593    
## subject126      0.296      1.879    0.16  0.87516    
## subject127     -0.172      1.286   -0.13  0.89357    
## subject128      2.413      1.354    1.78  0.07679 .  
## subject132     -3.123      1.631   -1.91  0.05733 .  
## subject133     -1.107      1.599   -0.69  0.48980    
## subject134     -2.431      1.817   -1.34  0.18306    
## subject135     -6.009      1.645   -3.65  0.00035 ***
## subject136     -0.311      1.802   -0.17  0.86333    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2.67 on 156 degrees of freedom
##   (6 observations deleted due to missingness)
## Multiple R-squared:  0.615,	Adjusted R-squared:  0.531 
## F-statistic: 7.32 on 34 and 156 DF,  p-value: <2e-16
```

From the first simple model, only men.stat is not significant. However, let's keep it for now while we test interactions. The next model is very complex with all 2-way, 3-way and 4-way interactions:

```r
fit.lab.2 <- lm(LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject + # single terms
                  tan.br.dr:age.sampling + tan.br.dr:men.stat + tan.br.dr:subject + # 2-way interactions
                  age.sampling:men.stat + age.sampling:subject + men.stat:subject + # 2-way interactions
                  tan.br.dr:age.sampling:men.stat + tan.br.dr:age.sampling:subject + # 3-way interactions
                  tan.br.dr:men.stat:subject + age.sampling:men.stat:subject + # 3-way interactions
                  tan.br.dr:age.sampling:men.stat:subject, # 4-way interaction 
                meta.girl.vag)
anova(fit.lab.2)
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##                                 Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr                        3    528   176.0   41.87 8.1e-13 ***
## age.sampling                     1     91    90.7   21.57 3.2e-05 ***
## men.stat                         1     26    25.9    6.17  0.0170 *  
## subject                         29   1132    39.0    9.28 8.7e-11 ***
## tan.br.dr:age.sampling           3     69    23.1    5.50  0.0027 ** 
## tan.br.dr:men.stat               3     26     8.6    2.04  0.1228    
## tan.br.dr:subject               41    357     8.7    2.07  0.0099 ** 
## age.sampling:men.stat            1      5     5.2    1.24  0.2708    
## age.sampling:subject            27    272    10.1    2.39  0.0051 ** 
## men.stat:subject                11     32     2.9    0.70  0.7340    
## tan.br.dr:age.sampling:men.stat  2     31    15.5    3.70  0.0330 *  
## tan.br.dr:age.sampling:subject  20    129     6.5    1.53  0.1186    
## age.sampling:men.stat:subject    5     10     2.1    0.50  0.7765    
## Residuals                       43    181     4.2                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.lab.2) # this is very long

anova(fit.lab.1, fit.lab.2) # compare 2 to 1
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:men.stat + 
##     age.sampling:subject + men.stat:subject + tan.br.dr:age.sampling:men.stat + 
##     tan.br.dr:age.sampling:subject + tan.br.dr:men.stat:subject + 
##     age.sampling:men.stat:subject + tan.br.dr:age.sampling:men.stat:subject
##   Res.Df  RSS  Df Sum of Sq    F Pr(>F)   
## 1    156 1113                             
## 2     43  181 113       932 1.96 0.0068 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Outcome of model 2 vs. 1: models are significantly different (p=0.0068), so keep more complex model 2 and reduce further in next step. Now we will drop the insignificant terms from fit.lab.2:

```r
fit.lab.3 <- lm(LAB.logit ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:age.sampling + tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject + # 2-way interactions
                  tan.br.dr:age.sampling:men.stat, # 3-way interaction
                meta.girl.vag)
anova(fit.lab.3)
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##                                 Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr                        3    528   176.0   38.25 2.3e-15 ***
## age.sampling                     1     91    90.7   19.71 2.9e-05 ***
## subject                         29   1117    38.5    8.37 2.9e-14 ***
## tan.br.dr:age.sampling           3     69    23.1    5.02  0.0031 ** 
## tan.br.dr:men.stat               4     67    16.7    3.63  0.0091 ** 
## tan.br.dr:subject               41    357     8.7    1.89  0.0077 ** 
## age.sampling:subject            27    277    10.2    2.23  0.0033 ** 
## tan.br.dr:age.sampling:men.stat  3     21     6.9    1.49  0.2233    
## Residuals                       79    364     4.6                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.lab.3) # this is very long
anova(fit.lab.2, fit.lab.3) # compare 3 to 2
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:men.stat + 
##     age.sampling:subject + men.stat:subject + tan.br.dr:age.sampling:men.stat + 
##     tan.br.dr:age.sampling:subject + tan.br.dr:men.stat:subject + 
##     age.sampling:men.stat:subject + tan.br.dr:age.sampling:men.stat:subject
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject + 
##     tan.br.dr:age.sampling:men.stat
##   Res.Df RSS  Df Sum of Sq    F Pr(>F)
## 1     43 181                          
## 2     79 364 -36      -183 1.21   0.28
```

Outcome of model 3 vs. 2: models are not significantly different (p=0.2754), so keep simpler model 3. Next, try dropping the 3-way term from fit.lab.3:

```r
fit.lab.4 <- lm(LAB.logit ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:age.sampling + tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject, # 2-way interactions
                meta.girl.vag)
anova(fit.lab.4)
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##                        Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr               3    528   176.0   37.58 2.2e-15 ***
## age.sampling            1     91    90.7   19.36 3.2e-05 ***
## subject                29   1117    38.5    8.22 2.5e-14 ***
## tan.br.dr:age.sampling  3     69    23.1    4.93  0.0034 ** 
## tan.br.dr:men.stat      4     67    16.7    3.56  0.0099 ** 
## tan.br.dr:subject      41    357     8.7    1.86  0.0088 ** 
## age.sampling:subject   27    277    10.2    2.19  0.0037 ** 
## Residuals              82    384     4.7                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.lab.4) # this is very long
anova(fit.lab.3, fit.lab.4) # compare 4 to 3
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject + 
##     tan.br.dr:age.sampling:men.stat
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject
##   Res.Df RSS Df Sum of Sq    F Pr(>F)
## 1     79 364                         
## 2     82 384 -3     -20.6 1.49   0.22
```

Outcome of model 4 vs. 3: models are not significantly different (p=0.2233), so keep simpler model 4. At this point all terms in model 4 are significant, but let's just try dropping the 2-way interaction with the highest p-value, tan.br.dr:men.stat :

```r
fit.lab.5 <- lm(LAB.logit ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:age.sampling + tan.br.dr:subject + age.sampling:subject, # 2-way interactions
                meta.girl.vag)
anova(fit.lab.5)
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##                        Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr               3    528   176.0   37.69 1.4e-15 ***
## age.sampling            1     91    90.7   19.41 3.1e-05 ***
## subject                29   1117    38.5    8.25 1.2e-14 ***
## tan.br.dr:age.sampling  3     69    23.1    4.94  0.0033 ** 
## tan.br.dr:subject      42    387     9.2    1.97  0.0041 ** 
## age.sampling:subject   27    300    11.1    2.38  0.0014 ** 
## Residuals              85    397     4.7                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.lab.5) # this is very long
anova(fit.lab.4, fit.lab.5) # compare 5 to 4
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:subject + age.sampling:subject
##   Res.Df RSS Df Sum of Sq    F Pr(>F)
## 1     82 384                         
## 2     85 397 -3     -12.9 0.92   0.44
```

Outcome of model 5 vs. 4: models are not significantly different (p=0.4358), so keep simpler model 5. Let's try dropping tan.br.dr:subject now:

```r
fit.lab.6 <- lm(LAB.logit ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:age.sampling + age.sampling:subject, # 2-way interactions
                meta.girl.vag)
anova(fit.lab.6)
```

```
## Analysis of Variance Table
## 
## Response: LAB.logit
##                         Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr                3    528   176.0   29.47 1.6e-14 ***
## age.sampling             1     91    90.7   15.18 0.00016 ***
## subject                 29   1117    38.5    6.45 4.6e-14 ***
## tan.br.dr:age.sampling   3     69    23.1    3.87 0.01098 *  
## age.sampling:subject    27    326    12.1    2.02 0.00502 ** 
## Residuals              127    759     6.0                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.lab.6) # this is fairly long
anova(fit.lab.5, fit.lab.6) # compare 6 to 4
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:subject + age.sampling:subject
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     age.sampling:subject
##   Res.Df RSS  Df Sum of Sq    F Pr(>F)   
## 1     85 397                             
## 2    127 759 -42      -362 1.84 0.0087 **
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Outcome of model 6 vs. 5: models are significantly different (p=0.0087), so we'll stick with slightly more complex model 5. Model 5 is considerably more complex than our simple model 1, but let's compare them to see how different they are:

```r
anova(fit.lab.1, fit.lab.5)
```

```
## Analysis of Variance Table
## 
## Model 1: LAB.logit ~ tan.br.dr + age.sampling + men.stat + subject
## Model 2: LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:subject + age.sampling:subject
##   Res.Df  RSS Df Sum of Sq    F  Pr(>F)    
## 1    156 1113                              
## 2     85  397 71       716 2.16 0.00037 ***
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Model 5 is significantly different from model 1 (p=0.0004). Final decision: keep model 5, which includes tan.br.dr, age.sampling, subject, and three 2-way interaction terms as predictors of LAB. View the summary of model 5, then output the ANOVA table:

```r
summary(fit.lab.5)
```

```
## 
## Call:
## lm(formula = LAB.logit ~ tan.br.dr + age.sampling + subject + 
##     tan.br.dr:age.sampling + tan.br.dr:subject + age.sampling:subject, 
##     data = meta.girl.vag)
## 
## Residuals:
##    Min     1Q Median     3Q    Max 
## -6.860 -0.639  0.000  0.623  5.383 
## 
## Coefficients: (47 not defined because of singularities)
##                           Estimate Std. Error t value Pr(>|t|)    
## (Intercept)                 9.5131    19.8590    0.48  0.63315    
## tan.br.dr.L                70.5174    20.9561    3.37  0.00115 ** 
## tan.br.dr.Q               -20.5858    15.1281   -1.36  0.17718    
## tan.br.dr.C               -24.0487    12.6263   -1.90  0.06021 .  
## age.sampling               -0.2299     1.5802   -0.15  0.88465    
## subject102                -12.6075    22.0536   -0.57  0.56905    
## subject103                -23.2484    33.7127   -0.69  0.49232    
## subject104                -29.3935    30.5100   -0.96  0.33808    
## subject105                -99.9794    28.8306   -3.47  0.00083 ***
## subject106                -35.7538    27.3902   -1.31  0.19530    
## subject107                -43.4793    26.8846   -1.62  0.10953    
## subject108                -44.5503    27.1920   -1.64  0.10504    
## subject109                -21.7842    39.5547   -0.55  0.58326    
## subject110                -10.1646     8.5815   -1.18  0.23953    
## subject111                -35.1729    25.0271   -1.41  0.16355    
## subject112                 31.9593    23.9361    1.34  0.18538    
## subject113                -10.4084    38.3996   -0.27  0.78701    
## subject114                -35.9722    28.8491   -1.25  0.21586    
## subject115                 18.4822    69.8050    0.26  0.79183    
## subject116               -115.3253   135.0250   -0.85  0.39545    
## subject118                -13.4977     2.5794   -5.23  1.2e-06 ***
## subject120                164.1614    71.2647    2.30  0.02369 *  
## subject121                -24.8671    45.9324   -0.54  0.58966    
## subject123                -63.0398    39.7439   -1.59  0.11642    
## subject124                -61.9129    30.5359   -2.03  0.04574 *  
## subject125               -188.0425   155.7046   -1.21  0.23052    
## subject126                  2.7053    60.9028    0.04  0.96467    
## subject127                 14.1043    53.0799    0.27  0.79110    
## subject128                -50.4998    42.4106   -1.19  0.23707    
## subject132                 82.3458    79.1062    1.04  0.30085    
## subject133               -178.6676    67.4875   -2.65  0.00967 ** 
## subject134                 19.8141    42.5784    0.47  0.64287    
## subject135                 -2.2796    53.0123   -0.04  0.96580    
## subject136                417.8299   144.9807    2.88  0.00500 ** 
## tan.br.dr.L:age.sampling   -5.6011     1.6461   -3.40  0.00102 ** 
## tan.br.dr.Q:age.sampling    1.8155     1.2255    1.48  0.14219    
## tan.br.dr.C:age.sampling    1.9479     1.0074    1.93  0.05648 .  
## tan.br.dr.L:subject102      8.3469     8.4204    0.99  0.32437    
## tan.br.dr.Q:subject102     -8.8192     6.3031   -1.40  0.16539    
## tan.br.dr.C:subject102          NA         NA      NA       NA    
## tan.br.dr.L:subject103     -2.9844     5.0268   -0.59  0.55429    
## tan.br.dr.Q:subject103     -2.1466     2.4649   -0.87  0.38629    
## tan.br.dr.C:subject103      2.4500     1.8191    1.35  0.18163    
## tan.br.dr.L:subject104     -6.7987     7.9577   -0.85  0.39531    
## tan.br.dr.Q:subject104      2.1540     5.9839    0.36  0.71977    
## tan.br.dr.C:subject104          NA         NA      NA       NA    
## tan.br.dr.L:subject105    -23.0047     8.6611   -2.66  0.00944 ** 
## tan.br.dr.Q:subject105          NA         NA      NA       NA    
## tan.br.dr.C:subject105          NA         NA      NA       NA    
## tan.br.dr.L:subject106      3.8661     6.8393    0.57  0.57338    
## tan.br.dr.Q:subject106     -5.9276     4.5282   -1.31  0.19405    
## tan.br.dr.C:subject106          NA         NA      NA       NA    
## tan.br.dr.L:subject107      6.8920     6.6222    1.04  0.30095    
## tan.br.dr.Q:subject107     -5.7931     4.8137   -1.20  0.23214    
## tan.br.dr.C:subject107          NA         NA      NA       NA    
## tan.br.dr.L:subject108      6.9704     6.0946    1.14  0.25596    
## tan.br.dr.Q:subject108      8.2850     5.3882    1.54  0.12786    
## tan.br.dr.C:subject108          NA         NA      NA       NA    
## tan.br.dr.L:subject109     -8.0313     6.3933   -1.26  0.21249    
## tan.br.dr.Q:subject109      0.8836     2.4619    0.36  0.72055    
## tan.br.dr.C:subject109      0.7943     2.1023    0.38  0.70651    
## tan.br.dr.L:subject110          NA         NA      NA       NA    
## tan.br.dr.Q:subject110          NA         NA      NA       NA    
## tan.br.dr.C:subject110          NA         NA      NA       NA    
## tan.br.dr.L:subject111     -3.0107     4.8411   -0.62  0.53567    
## tan.br.dr.Q:subject111      0.0696     2.6264    0.03  0.97892    
## tan.br.dr.C:subject111      0.2585     2.0428    0.13  0.89961    
## tan.br.dr.L:subject112    -24.8112     8.3122   -2.98  0.00370 ** 
## tan.br.dr.Q:subject112    -11.3428     4.9424   -2.30  0.02420 *  
## tan.br.dr.C:subject112          NA         NA      NA       NA    
## tan.br.dr.L:subject113     -3.9607     6.8811   -0.58  0.56641    
## tan.br.dr.Q:subject113      2.2848     5.9045    0.39  0.69976    
## tan.br.dr.C:subject113          NA         NA      NA       NA    
## tan.br.dr.L:subject114    -18.0065    10.3198   -1.74  0.08463 .  
## tan.br.dr.Q:subject114      7.2455     6.0571    1.20  0.23494    
## tan.br.dr.C:subject114          NA         NA      NA       NA    
## tan.br.dr.L:subject115    -27.6026    10.4904   -2.63  0.01010 *  
## tan.br.dr.Q:subject115          NA         NA      NA       NA    
## tan.br.dr.C:subject115          NA         NA      NA       NA    
## tan.br.dr.L:subject116     -5.6843     7.3692   -0.77  0.44263    
## tan.br.dr.Q:subject116          NA         NA      NA       NA    
## tan.br.dr.C:subject116          NA         NA      NA       NA    
## tan.br.dr.L:subject118          NA         NA      NA       NA    
## tan.br.dr.Q:subject118          NA         NA      NA       NA    
## tan.br.dr.C:subject118          NA         NA      NA       NA    
## tan.br.dr.L:subject120      9.0691     9.0362    1.00  0.31840    
## tan.br.dr.Q:subject120      1.3282     5.4807    0.24  0.80910    
## tan.br.dr.C:subject120          NA         NA      NA       NA    
## tan.br.dr.L:subject121    -25.6551     8.7291   -2.94  0.00424 ** 
## tan.br.dr.Q:subject121    -12.2855     7.3042   -1.68  0.09624 .  
## tan.br.dr.C:subject121          NA         NA      NA       NA    
## tan.br.dr.L:subject123    -15.5585     5.9598   -2.61  0.01068 *  
## tan.br.dr.Q:subject123          NA         NA      NA       NA    
## tan.br.dr.C:subject123          NA         NA      NA       NA    
## tan.br.dr.L:subject124    -24.0409     8.5608   -2.81  0.00618 ** 
## tan.br.dr.Q:subject124     -7.7261     5.7934   -1.33  0.18590    
## tan.br.dr.C:subject124          NA         NA      NA       NA    
## tan.br.dr.L:subject125    -21.7414    13.5635   -1.60  0.11266    
## tan.br.dr.Q:subject125          NA         NA      NA       NA    
## tan.br.dr.C:subject125          NA         NA      NA       NA    
## tan.br.dr.L:subject126     -5.7340    13.5347   -0.42  0.67289    
## tan.br.dr.Q:subject126          NA         NA      NA       NA    
## tan.br.dr.C:subject126          NA         NA      NA       NA    
## tan.br.dr.L:subject127      6.4943     8.8991    0.73  0.46754    
## tan.br.dr.Q:subject127          NA         NA      NA       NA    
## tan.br.dr.C:subject127          NA         NA      NA       NA    
## tan.br.dr.L:subject128    -15.0513     9.7873   -1.54  0.12780    
## tan.br.dr.Q:subject128          NA         NA      NA       NA    
## tan.br.dr.C:subject128          NA         NA      NA       NA    
## tan.br.dr.L:subject132     18.4158    12.9038    1.43  0.15719    
## tan.br.dr.Q:subject132          NA         NA      NA       NA    
## tan.br.dr.C:subject132          NA         NA      NA       NA    
## tan.br.dr.L:subject133     -7.4674     6.2519   -1.19  0.23563    
## tan.br.dr.Q:subject133          NA         NA      NA       NA    
## tan.br.dr.C:subject133          NA         NA      NA       NA    
## tan.br.dr.L:subject134          NA         NA      NA       NA    
## tan.br.dr.Q:subject134          NA         NA      NA       NA    
## tan.br.dr.C:subject134          NA         NA      NA       NA    
## tan.br.dr.L:subject135          NA         NA      NA       NA    
## tan.br.dr.Q:subject135          NA         NA      NA       NA    
## tan.br.dr.C:subject135          NA         NA      NA       NA    
## tan.br.dr.L:subject136     31.6250    14.0760    2.25  0.02725 *  
## tan.br.dr.Q:subject136          NA         NA      NA       NA    
## tan.br.dr.C:subject136          NA         NA      NA       NA    
## age.sampling:subject102     0.4601     1.8210    0.25  0.80112    
## age.sampling:subject103     1.7992     2.7345    0.66  0.51234    
## age.sampling:subject104     2.2294     2.4587    0.91  0.36709    
## age.sampling:subject105     8.4202     2.4927    3.38  0.00110 ** 
## age.sampling:subject106     2.7289     2.1656    1.26  0.21107    
## age.sampling:subject107     3.1289     2.0954    1.49  0.13908    
## age.sampling:subject108     3.9080     2.1646    1.81  0.07455 .  
## age.sampling:subject109     1.9351     3.3649    0.58  0.56676    
## age.sampling:subject110         NA         NA      NA       NA    
## age.sampling:subject111     3.0282     2.0797    1.46  0.14905    
## age.sampling:subject112    -3.8767     1.9887   -1.95  0.05454 .  
## age.sampling:subject113     0.8088     3.0682    0.26  0.79273    
## age.sampling:subject114     3.5293     2.4985    1.41  0.16143    
## age.sampling:subject115    -2.1423     5.9227   -0.36  0.71847    
## age.sampling:subject116     9.4608    11.5523    0.82  0.41511    
## age.sampling:subject118         NA         NA      NA       NA    
## age.sampling:subject120   -13.6958     5.7471   -2.38  0.01940 *  
## age.sampling:subject121     1.3883     4.0943    0.34  0.73538    
## age.sampling:subject123     5.2799     3.4811    1.52  0.13304    
## age.sampling:subject124     5.1021     2.6555    1.92  0.05805 .  
## age.sampling:subject125    16.3032    13.4737    1.21  0.22963    
## age.sampling:subject126    -0.6217     5.4771   -0.11  0.90990    
## age.sampling:subject127    -1.0808     4.1790   -0.26  0.79655    
## age.sampling:subject128     4.3808     3.6084    1.21  0.22809    
## age.sampling:subject132    -6.8804     6.6488   -1.03  0.30369    
## age.sampling:subject133    14.5895     5.5981    2.61  0.01081 *  
## age.sampling:subject134    -2.0977     3.7789   -0.56  0.58027    
## age.sampling:subject135    -0.5324     4.5812   -0.12  0.90775    
## age.sampling:subject136   -37.1749    12.8069   -2.90  0.00471 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 2.16 on 85 degrees of freedom
##   (6 observations deleted due to missingness)
## Multiple R-squared:  0.863,	Adjusted R-squared:  0.693 
## F-statistic: 5.08 on 105 and 85 DF,  p-value: 1.5e-13
```

```r
write.csv(anova(fit.lab.5), file="anova-lab-best-lm.csv")
```

## Model comparisons for vaginal pH

Now we use the same approach to find a multiple linear regression model for vaginal pH. First, a very simple model with four variables and no interactions:

```r
fit.ph.1 <- lm(ph ~ tan.br.dr + age.sampling + men.stat + subject, meta.girl.vag)
anova(fit.ph.1) 
```

```
## Analysis of Variance Table
## 
## Response: ph
##              Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr     3   16.6    5.54   13.36 2.2e-07 ***
## age.sampling  1   15.8   15.84   38.23 1.4e-08 ***
## men.stat      1    0.1    0.07    0.17    0.68    
## subject      24   75.8    3.16    7.62 1.2e-13 ***
## Residuals    99   41.0    0.41                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(fit.ph.1)
```

```
## 
## Call:
## lm(formula = ph ~ tan.br.dr + age.sampling + men.stat + subject, 
##     data = meta.girl.vag)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.6817 -0.2764 -0.0298  0.2292  2.3463 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    7.6945     1.8212    4.22  5.3e-05 ***
## tan.br.dr.L   -0.3665     0.2478   -1.48  0.14232    
## tan.br.dr.Q    0.2340     0.1477    1.58  0.11639    
## tan.br.dr.C   -0.0516     0.1251   -0.41  0.68095    
## age.sampling  -0.2817     0.1426   -1.98  0.05094 .  
## men.statpost  -0.1784     0.2167   -0.82  0.41230    
## subject102     2.0139     0.4191    4.81  5.5e-06 ***
## subject103     1.6437     0.4077    4.03  0.00011 ***
## subject104     2.5991     0.5079    5.12  1.5e-06 ***
## subject105     2.0805     0.4966    4.19  6.1e-05 ***
## subject106     0.7946     0.4411    1.80  0.07469 .  
## subject107     0.8045     0.4104    1.96  0.05278 .  
## subject108     0.5608     0.4210    1.33  0.18593    
## subject109     0.9266     0.4370    2.12  0.03646 *  
## subject111     0.7551     0.4625    1.63  0.10574    
## subject112     0.0469     0.4849    0.10  0.92321    
## subject113     0.3417     0.4713    0.72  0.47019    
## subject114     0.3815     0.5249    0.73  0.46913    
## subject116     2.1704     0.7626    2.85  0.00538 ** 
## subject120     2.1808     0.4885    4.46  2.1e-05 ***
## subject121     0.8778     0.5243    1.67  0.09722 .  
## subject123     0.7250     0.6337    1.14  0.25533    
## subject124     0.0648     0.5909    0.11  0.91288    
## subject127     0.0725     0.7388    0.10  0.92205    
## subject128     0.5120     0.5101    1.00  0.31796    
## subject132     1.9973     0.5145    3.88  0.00019 ***
## subject133     0.1714     0.5191    0.33  0.74191    
## subject134    -0.0423     0.7609   -0.06  0.95579    
## subject135     1.1929     0.6097    1.96  0.05320 .  
## subject136     1.5061     0.5277    2.85  0.00525 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.644 on 99 degrees of freedom
##   (68 observations deleted due to missingness)
## Multiple R-squared:  0.725,	Adjusted R-squared:  0.645 
## F-statistic: 9.01 on 29 and 99 DF,  p-value: <2e-16
```

Once again, only men.stat is not significant. However, let's keep it for now while we test interactions. The next is a complex model with all 2-way, 3-way and 4-way interactions:

```r
fit.ph.2 <- lm(ph ~ tan.br.dr + age.sampling + men.stat + subject + # single terms
                  tan.br.dr:age.sampling + tan.br.dr:men.stat + tan.br.dr:subject + # 2-way interactions
                  age.sampling:men.stat + age.sampling:subject + men.stat:subject + # 2-way interactions
                  tan.br.dr:age.sampling:men.stat + tan.br.dr:age.sampling:subject + # 3-way interactions
                  tan.br.dr:men.stat:subject + age.sampling:men.stat:subject + # 3-way interactions
                  tan.br.dr:age.sampling:men.stat:subject, # 4-way interaction 
                meta.girl.vag)
anova(fit.ph.2)
```

```
## Analysis of Variance Table
## 
## Response: ph
##                                 Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr                        3   16.6    5.54   25.67 1.2e-07 ***
## age.sampling                     1   15.8   15.84   73.49 9.1e-09 ***
## men.stat                         1    0.1    0.07    0.33  0.5719    
## subject                         24   75.8    3.16   14.64 3.2e-09 ***
## tan.br.dr:age.sampling           3    0.7    0.23    1.06  0.3842    
## tan.br.dr:men.stat               3    2.3    0.77    3.56  0.0290 *  
## tan.br.dr:subject               27   13.5    0.50    2.32  0.0202 *  
## age.sampling:men.stat            1    0.6    0.56    2.58  0.1215    
## age.sampling:subject            20   10.3    0.52    2.40  0.0214 *  
## men.stat:subject                 7    5.4    0.77    3.57  0.0089 ** 
## tan.br.dr:age.sampling:men.stat  1    0.3    0.25    1.16  0.2918    
## tan.br.dr:age.sampling:subject  11    2.7    0.25    1.15  0.3721    
## age.sampling:men.stat:subject    2    0.1    0.05    0.21  0.8128    
## Residuals                       24    5.2    0.22                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.ph.2) # this is very long
anova(fit.ph.1, fit.ph.2) # compare 2 to 1
```

```
## Analysis of Variance Table
## 
## Model 1: ph ~ tan.br.dr + age.sampling + men.stat + subject
## Model 2: ph ~ tan.br.dr + age.sampling + men.stat + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:men.stat + 
##     age.sampling:subject + men.stat:subject + tan.br.dr:age.sampling:men.stat + 
##     tan.br.dr:age.sampling:subject + tan.br.dr:men.stat:subject + 
##     age.sampling:men.stat:subject + tan.br.dr:age.sampling:men.stat:subject
##   Res.Df  RSS Df Sum of Sq    F Pr(>F)  
## 1     99 41.0                           
## 2     24  5.2 75      35.9 2.22  0.015 *
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Outcome of model 2 vs. 1: models are significantly different (p=0.0154), so keep more complex model 2 and reduce further. Next, drop the insignificant terms from fit.ph.2:

```r
fit.ph.3 <- lm(ph ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject + men.stat:subject, # 2-way interactions 
                meta.girl.vag)
anova(fit.ph.3)
```

```
## Analysis of Variance Table
## 
## Response: ph
##                      Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr             3   16.6    5.54   18.46 8.6e-08 ***
## age.sampling          1   15.8   15.84   52.84 6.0e-09 ***
## subject              24   75.6    3.15   10.50 4.6e-11 ***
## tan.br.dr:men.stat    4    3.0    0.74    2.48   0.059 .  
## tan.br.dr:subject    27   13.3    0.49    1.64   0.073 .  
## age.sampling:subject 20    9.1    0.46    1.53   0.123    
## subject:men.stat      7    3.3    0.47    1.58   0.169    
## Residuals            42   12.6    0.30                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.ph.3) # this is very long
anova(fit.ph.2, fit.ph.3) # compare 3 to 2
```

```
## Analysis of Variance Table
## 
## Model 1: ph ~ tan.br.dr + age.sampling + men.stat + subject + tan.br.dr:age.sampling + 
##     tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:men.stat + 
##     age.sampling:subject + men.stat:subject + tan.br.dr:age.sampling:men.stat + 
##     tan.br.dr:age.sampling:subject + tan.br.dr:men.stat:subject + 
##     age.sampling:men.stat:subject + tan.br.dr:age.sampling:men.stat:subject
## Model 2: ph ~ tan.br.dr + age.sampling + subject + tan.br.dr:men.stat + 
##     tan.br.dr:subject + age.sampling:subject + men.stat:subject
##   Res.Df   RSS  Df Sum of Sq    F Pr(>F)  
## 1     24  5.17                            
## 2     42 12.59 -18     -7.42 1.91  0.069 .
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

Outcome of model 3 vs. 2: models 2 and 3 are not significantly different (p=0.0691), so keep simpler model 3. Next, drop the insignificant terms from fit.ph.3:

```r
fit.ph.4 <- lm(ph ~ tan.br.dr + age.sampling + subject + # single terms
                  tan.br.dr:men.stat + tan.br.dr:subject, # 2-way interactions 
                meta.girl.vag)
anova(fit.ph.4)
```

```
## Analysis of Variance Table
## 
## Response: ph
##                    Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr           3   16.6    5.54   15.25 1.0e-07 ***
## age.sampling        1   15.8   15.84   43.64 6.8e-09 ***
## subject            24   75.6    3.15    8.67 8.2e-13 ***
## tan.br.dr:men.stat  4    3.0    0.74    2.04   0.098 .  
## tan.br.dr:subject  27   13.3    0.49    1.36   0.156    
## Residuals          69   25.1    0.36                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
# summary(fit.ph.4) # this is fairly long
anova(fit.ph.3, fit.ph.4) # compare 4 to 3
```

```
## Analysis of Variance Table
## 
## Model 1: ph ~ tan.br.dr + age.sampling + subject + tan.br.dr:men.stat + 
##     tan.br.dr:subject + age.sampling:subject + men.stat:subject
## Model 2: ph ~ tan.br.dr + age.sampling + subject + tan.br.dr:men.stat + 
##     tan.br.dr:subject
##   Res.Df  RSS  Df Sum of Sq    F Pr(>F)
## 1     42 12.6                          
## 2     69 25.1 -27     -12.5 1.54    0.1
```

Outcome of model 4 vs. 3: models 3 and 4 are not significantly different (p=0.1026), so keep simpler model 4. Next, drop the insignificant terms from fit.ph.4 (leaving three single terms only):

```r
fit.ph.5 <- lm(ph ~ tan.br.dr + age.sampling + subject, meta.girl.vag)
anova(fit.ph.5)
```

```
## Analysis of Variance Table
## 
## Response: ph
##               Df Sum Sq Mean Sq F value  Pr(>F)    
## tan.br.dr      3   16.6    5.54   13.40 2.0e-07 ***
## age.sampling   1   15.8   15.84   38.35 1.3e-08 ***
## subject       24   75.6    3.15    7.62 1.0e-13 ***
## Residuals    100   41.3    0.41                    
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
```

```r
summary(fit.ph.5)
```

```
## 
## Call:
## lm(formula = ph ~ tan.br.dr + age.sampling + subject, data = meta.girl.vag)
## 
## Residuals:
##     Min      1Q  Median      3Q     Max 
## -1.7757 -0.2916 -0.0326  0.2354  2.3551 
## 
## Coefficients:
##              Estimate Std. Error t value Pr(>|t|)    
## (Intercept)    8.4647     1.5601    5.43  4.0e-07 ***
## tan.br.dr.L   -0.3632     0.2474   -1.47  0.14511    
## tan.br.dr.Q    0.2207     0.1466    1.51  0.13534    
## tan.br.dr.C   -0.0376     0.1238   -0.30  0.76177    
## age.sampling  -0.3449     0.1199   -2.88  0.00492 ** 
## subject102     1.8800     0.3856    4.88  4.1e-06 ***
## subject103     1.5879     0.4013    3.96  0.00014 ***
## subject104     2.4001     0.4460    5.38  4.9e-07 ***
## subject105     2.0031     0.4868    4.12  8.0e-05 ***
## subject106     0.7139     0.4294    1.66  0.09952 .  
## subject107     0.7558     0.4055    1.86  0.06526 .  
## subject108     0.5570     0.4203    1.33  0.18811    
## subject109     0.8335     0.4214    1.98  0.05069 .  
## subject111     0.6544     0.4453    1.47  0.14488    
## subject112    -0.1354     0.4308   -0.31  0.75401    
## subject113     0.1915     0.4339    0.44  0.65990    
## subject114     0.2246     0.4884    0.46  0.64653    
## subject116     2.1391     0.7604    2.81  0.00591 ** 
## subject120     2.0727     0.4698    4.41  2.6e-05 ***
## subject121     0.8052     0.5160    1.56  0.12177    
## subject123     0.6293     0.6219    1.01  0.31401    
## subject124    -0.1915     0.5014   -0.38  0.70327    
## subject127     0.1129     0.7360    0.15  0.87844    
## subject128     0.4128     0.4949    0.83  0.40620    
## subject132     1.9235     0.5058    3.80  0.00025 ***
## subject133     0.1786     0.5181    0.34  0.73104    
## subject134    -0.0540     0.7595   -0.07  0.94351    
## subject135     1.1892     0.6087    1.95  0.05353 .  
## subject136     1.4370     0.5201    2.76  0.00682 ** 
## ---
## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
## 
## Residual standard error: 0.643 on 100 degrees of freedom
##   (68 observations deleted due to missingness)
## Multiple R-squared:  0.723,	Adjusted R-squared:  0.646 
## F-statistic: 9.34 on 28 and 100 DF,  p-value: <2e-16
```

```r
anova(fit.ph.4, fit.ph.5) # compare 5 to 4
```

```
## Analysis of Variance Table
## 
## Model 1: ph ~ tan.br.dr + age.sampling + subject + tan.br.dr:men.stat + 
##     tan.br.dr:subject
## Model 2: ph ~ tan.br.dr + age.sampling + subject
##   Res.Df  RSS  Df Sum of Sq    F Pr(>F)
## 1     69 25.1                          
## 2    100 41.3 -31     -16.3 1.44    0.1
```

Outcome of model 5 vs. 4: models 4 and 5 are not significantly different (p=0.1036), so keep simpler model 5. All three factors in model 5 are highly significant, so it doesn't make sense to try taking another term out. Finally, compare model 5 to the first simple model 1:

```r
anova(fit.ph.1, fit.ph.5)
```

```
## Analysis of Variance Table
## 
## Model 1: ph ~ tan.br.dr + age.sampling + men.stat + subject
## Model 2: ph ~ tan.br.dr + age.sampling + subject
##   Res.Df  RSS Df Sum of Sq    F Pr(>F)
## 1     99 41.0                         
## 2    100 41.3 -1    -0.281 0.68   0.41
```

Final decision: models 1 and 5 are not significantly different (p=0.4123), so keep simpler model 5. Output the ANOVA table:

```r
write.csv(anova(fit.ph.5), file="table-s7-anova-ph-best-lm.csv")
```

## Recap of regression model comparisons
To recap, we selected the following as optimal multiple linear regression models for LAB (logit-transformed) and vaginal pH:

* LAB.logit ~ tan.br.dr + age.sampling + subject + tan.br.dr:age.sampling + tan.br.dr:men.stat + tan.br.dr:subject + age.sampling:subject
* ph ~ tan.br.dr + age.sampling + subject

***
# Part III. Plotting trends in lactic acid bacteria and vaginal pH

## Trends in LAB associated with participant metadata

```r
## Boxplot of LAB.logit ~ Tanner breast (significant in linear model)
gg.lab.logit.tb <- ggplot(subset(meta.girl.vag, tan.br.dr %in% c(2,3,4,5)), 
                          aes(y=LAB.logit, x=tan.br.dr, fill=tan.br.dr))
box.lab.logit.tb <- gg.lab.logit.tb + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.tanner[2:5], name="Tanner\nBreast") +
  guides(fill=FALSE) +
  xlab("Tanner breast") +
  ylab("Logit-transformed LAB proportion") +
  theme_cust

print(box.lab.logit.tb)
```

![plot of chunk lab-trends](./adolescent-supp-03_files/figure-html/lab-trends1.png) 

```r
## Boxplot of LAB.logit ~ menarche status (not significant in linear model)
gg.lab.logit.ms <- ggplot(subset(meta.girl.vag, men.stat %in% c("pre","post")), 
                          aes(y=LAB.logit, x=men.stat, fill=men.stat))
box.lab.logit.ms <- gg.lab.logit.ms + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.men.stat[c(3,2)], name="Menarche\nStatus") +
  guides(fill=FALSE) +
  xlab("Menarche status") +
  ylab("Logit-transformed LAB proportion") +
  theme_cust

print(box.lab.logit.ms)
```

![plot of chunk lab-trends](./adolescent-supp-03_files/figure-html/lab-trends2.png) 

```r
## Scatterplot of LAB.logit ~ age (significant in linear model)
gg.lab.logit.age <- ggplot(meta.girl.vag, 
                           aes(y=LAB.logit, x=age.sampling))
scatter.lab.logit.age <- gg.lab.logit.age + 
  geom_point(size=2, alpha=0.7) +
  xlab("Age (y)") +
  ylab("Logit-transformed LAB proportion") +
  theme_cust

print(scatter.lab.logit.age)
```

![plot of chunk lab-trends](./adolescent-supp-03_files/figure-html/lab-trends3.png) 

```r
## Boxplot of LAB.logit ~ subject (significant in linear model)
gg.lab.logit.sub <- ggplot(meta.girl.vag, 
                           aes(y=LAB.logit, x=subject, fill=subject))
box.lab.logit.sub <- gg.lab.logit.sub + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.gm.pair[1:31], name="Subject") +
  guides(fill=FALSE) +
  xlab("Subject") +
  ylab("Logit-transformed LAB proportion") +
  theme_cust

print(box.lab.logit.sub)
```

![plot of chunk lab-trends](./adolescent-supp-03_files/figure-html/lab-trends4.png) 

## Trends in vaginal pH associated with participant metadata

```r
## Boxplot of pH ~ Tanner breast (significant in linear model)
gg.ph.tb <- ggplot(subset(meta.girl.vag, tan.br.dr %in% c(2,3,4,5)), 
                   aes(y=ph, x=tan.br.dr, fill=tan.br.dr))
box.ph.tb <- gg.ph.tb + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.tanner[2:5], name="Tanner\nBreast") +
  guides(fill=FALSE) +
  xlab("Tanner breast") +
  ylab("Vaginal pH") +
  theme_cust

print(box.ph.tb)
```

```
## Warning: Removed 62 rows containing non-finite values (stat_boxplot).
```

![plot of chunk ph-trends](./adolescent-supp-03_files/figure-html/ph-trends1.png) 

```r
## Boxplot of pH ~ menarche status (not significant in linear model)
gg.ph.ms <- ggplot(subset(meta.girl.vag, men.stat %in% c("pre","post")), 
                   aes(y=ph, x=men.stat, fill=men.stat))
box.ph.ms <- gg.ph.ms + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.men.stat[c(3,2)], name="Menarche\nStatus") +
  guides(fill=FALSE) +
  xlab("Menarche status") +
  ylab("Vaginal pH") +
  theme_cust

print(box.ph.ms)
```

```
## Warning: Removed 66 rows containing non-finite values (stat_boxplot).
```

![plot of chunk ph-trends](./adolescent-supp-03_files/figure-html/ph-trends2.png) 

```r
## Scatterplot of pH ~ age (significant in linear model)
gg.ph.age <- ggplot(meta.girl.vag, 
                    aes(y=ph, x=age.sampling))
scatter.ph.age <- gg.ph.age + 
  geom_point(size=2, alpha=0.7) +
  xlab("Age (y)") +
  ylab("Vaginal pH") +
  theme_cust

print(scatter.ph.age)
```

```
## Warning: Removed 67 rows containing missing values (geom_point).
```

![plot of chunk ph-trends](./adolescent-supp-03_files/figure-html/ph-trends3.png) 

```r
## Boxplot of ph ~ subject (significant in linear model)
gg.ph.sub <- ggplot(meta.girl.vag, 
                    aes(y=ph, x=subject, fill=subject))
box.ph.sub <- gg.ph.sub + 
  geom_boxplot(width=0.5, outlier.shape=1) +
  scale_fill_manual(values=col.gm.pair[1:31], name="Subject") +
  guides(fill=FALSE) +
  xlab("Subject") +
  ylab("Vaginal pH") +
  theme_cust

print(box.ph.sub)
```

```
## Warning: Removed 67 rows containing non-finite values (stat_boxplot).
```

![plot of chunk ph-trends](./adolescent-supp-03_files/figure-html/ph-trends4.png) 

## Figure 6. Trends in relative abundance of lactic acid bacteria and vaginal pH in relation to pubertal development and menarche status.
198 vaginal swabs were collected from 31 girls over time. Upper and lower panels show box plots of (a) the logit-transformed proportion of lactic acid bacteria (LAB; includes _Lactobacillus_, _Streptococcus_, _Aerococcus_ and _Facklamia_) and (b) vaginal pH. In the left column, box plots show the relationship to Tanner breast score (191 data points for LAB and 130 data points for pH out of 198 total samples). Box plots in the right column show the relationship to menarche status (197 data points for LAB and 130 data points for pH). In each plot the rectangular box represents the interquartile range, the whiskers represent the upper and lower quartiles, the horizontal line represents the median, and open circles represent outliers.


```r
## Make compound figure for manuscript with LAB.logit panels on top, pH on bottom
fig.lab.logit.ph.1 <- box.lab.logit.tb +
  ggtitle("a") +
  theme(plot.title=element_text(size=22, hjust=-0.15, vjust=1.5))
fig.lab.logit.ph.2 <- box.lab.logit.ms +
  ggtitle("") +
  ylab("") +
  theme(plot.title=element_text(size=22, hjust=-0.15, vjust=1.5))
fig.lab.logit.ph.3 <- box.ph.tb +
  ggtitle("b") +
  theme(plot.title=element_text(size=22, hjust=-0.15, vjust=1.5))
fig.lab.logit.ph.4 <- box.ph.ms +
  ggtitle("") +
  ylab("") +
  theme(plot.title=element_text(size=22, hjust=-0.15, vjust=1.5))

# pdf("fig-6-lab-logit-ph-trends.pdf", width=6, height=8)
multiplot(fig.lab.logit.ph.1, fig.lab.logit.ph.2,  
          fig.lab.logit.ph.3, fig.lab.logit.ph.4,  
          layout=matrix(c(1,1,1,2,2,3,3,3,4,4), ncol=5, byrow=TRUE))
```

```
## Warning: Removed 62 rows containing non-finite values (stat_boxplot).
## Warning: Removed 66 rows containing non-finite values (stat_boxplot).
```

![plot of chunk fig-6-lab-ph-trends](./adolescent-supp-03_files/figure-html/fig-6-lab-ph-trends.png) 

```r
dev.off()
```

```
## null device 
##           1
```

## Discordance between high LAB and low pH

An interesting observation was that not all samples with high proportions of LAB were associated with a low pH. We hypothesized this could be due to lower total bacterial counts in the vaginas of adolescent girls. We performed a coarse test of that hypothesis (see Supplemental File 5, 'adolescent-supp-05.Rmd) but were unable to detect a significant difference in estimated number of 16S rRNA gene copies.

## Figure 7. Relationship between proportion of lactic acid bacteria and pH in vaginal samples collected from girls.
197 vaginal microbiota samples from girls are plotted to show the relationship between the proportion of lactic acid bacteria (LAB; includes _Lactobacillus_, _Streptococcus_, _Aerococcus_ and _Facklamia_) on the x-axis and vaginal pH on the y-axis. Menarche status is either premenarcheal (open circles) or postmenarcheal (filled circles). Points are color-coded to indicate Tanner breast stage as indicated by the legend at right (NA values are colored gray). Points are slightly jittered to decrease crowding around similar values, but not so much as to distort interpretation of the data. The magenta dashed line at vaginal pH 4.5 represents the upper range of the traditional ‘hallmark’ healthy vaginal pH of 4.0-4.5.


```r
## LAB vs. ph
gg.lab.ph <- ggplot(subset(meta, type=="girl" & site=="vag"), 
                    aes(x=LAB, y=ph, col=factor(tan.br.dr), shape=men.stat))

gg.lab.ph + geom_hline(yintercept=4.5, linetype=2, color="#8B1C62", alpha=0.8) +
  geom_jitter(position=position_jitter(width = 0.01, height = 0.05), size=2.5, alpha=0.8) +
  scale_color_manual(values=col.tanner, name="Tanner\nBreast", na.value="gray70") +
  scale_shape_manual(values=c(16,1), name="Menarche\nStatus", breaks=c("pre", "post")) +
  xlab("Proportion of LAB") +
  ylab("Vaginal pH") +
  theme_cust
```

```
## Warning: Removed 68 rows containing missing values (geom_point).
```

![plot of chunk fig-7-lab-ph](./adolescent-supp-03_files/figure-html/fig-7-lab-ph.png) 

```r
# Uncomment line below to save as PDF
# ggsave("fig-7-lab-ph-scatterplot.pdf", width=8, height=6, units="in")
```

End of Supplemental File 3. For a description of qPCR performed to test whether low-pH samples had more bacterial 16S gene copies than high-pH samples, see *Supplemental File 5. Comparison of low-pH vs. high-pH vaginal microbiota with pan-bacterial 16S rRNA qPCR assay.* For more exploration of the vulvar microbiota, see *Supplemental File 6. Comparisons of vaginal and vulvar microbiota of girls.*
