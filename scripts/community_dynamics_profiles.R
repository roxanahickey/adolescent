load("data-postproc/02-hclust-pcoa-last-run.RData")

library(ggplot2)
library(reshape)
library(gridExtra)

## g_legend function from http://stackoverflow.com/questions/11883844/
## inserting-a-table-under-the-legend-in-a-ggplot2-histogram
g_legend <- function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

## Appendix 2: Individual summaries of vaginal and vulvar microbiota composition and associated metadata for all study participants

pdf(file="supplemental/appendix-s2.pdf", 
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
dev.off()