##########################################################################
library(reshape)
library(ggplot2)
library(scales)
library(grid)
library(gridExtra)

setwd("~/Documents/research/pg-adolescent/analysis-rjh/git-adolescent/")
load("data-postproc/01-data-prep-last-run.RData")

g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

pdf(file="summary-plots-by-participant.pdf", 
    width=8.5, height=11, pointsize=8)
for (i in as.numeric(unique(meta$gm.pair))){
  pick <- which(as.numeric(meta$gm.pair) == i)
  if (length(pick) > 1){
    prop.tmp <- prop.red[,pick]
    meta.tmp <- meta[pick,]
    
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
                
    ## Community composition
    gg.taxa <- ggplot(meta.prop25.lg.tmp, 
                     aes(x=visit, y=prop, fill=taxon))
    gg.taxa.bar <- gg.taxa + geom_bar(stat="identity") +
      facet_wrap(~ type.site, ncol=1, scale="fixed") +
      scale_fill_manual(values=col.taxa[taxa.pick25], name="Taxon") +
      scale_x_discrete(limits=seq(1, max(meta.prop25.lg.tmp$visit), by=1), 
                       name="") +
      ylab("Proportion") +
      ggtitle(paste(unique(meta.tmp$gm.pair.fullID))) +
      theme_cust +
      theme(axis.text.x=element_blank())
            
    ## Now we need to go back to the original metadata table which has each
    ## visit listed in its own row. This way we don't have duplicate menarche
    ## status and Tanner scores to plot (one matching vagina, one matching
    ## vulva when they are identical)
    meta2.tmp <- subset(meta.orig, 
                        sample.ID.vag %in% rownames(meta.tmp) | 
                          sample.ID.vul %in% rownames(meta.tmp))

    ## Menarche status
    gg.men.stat <- ggplot(subset(meta2.tmp, type=="girl"), 
                          aes(x=visit, y=subject, color=men.stat))
    gg.men.stat.dot <- gg.men.stat + 
      geom_point(size=7) +  
      scale_color_manual(values=col.men.stat[2:3],
                         breaks=c("pre","post"),
                         na.value="gray70",
                         name="Menarche\nStatus") +
      scale_x_discrete(limits=seq(1, max(meta.prop25.lg.tmp$visit), by=1), 
                       name="") +
      scale_y_discrete(labels="Men", name="") +
      theme_cust +
      theme(axis.text.x=element_blank())
        
    ## Tanner
    tanner.tmp <- meta2.tmp[meta2.tmp$type=="girl",
                            c("visit","tan.br.dr","tan.gen.dr")]
    tanner.lg.tmp <- melt(tanner.tmp, id.vars="visit")
    
    gg.tanner <- ggplot(tanner.lg.tmp, 
                        aes(x=visit, y=variable, color=factor(value)))
    gg.tanner.dot <- gg.tanner +  
      geom_point(size=7) +  
      scale_color_manual(values=col.tanner, 
                         name="Tanner\nScore", 
                         na.value="gray70") +  
      scale_x_discrete(1:14, name="Visit") +
      scale_y_discrete(breaks=c("tan.gen.dr","tan.br.dr"), 
                       labels=c("Gen", "Bre"), name="") +
      theme_cust
    
#     # pH
#     gg.ph <- ggplot(subset(meta2.tmp, type=="girl"), 
#                     aes(x=visit, y=ph))
#     gg.ph.line <- gg.ph +
#       geom_point(size=2) +
#       geom_line() +
#       scale_x_discrete(1:max(meta.prop25.lg.tmp$visit), 
#                        name="") +
#       scale_y_continuous(limits=c(4.0, 7.5)) +
#       geom_text(label=meta2.tmp$ph[meta2.tmp$type=="girl"], 
#                 vjust=-0.5) +
#       ylab("Vaginal pH") +
#       theme_cust_nominor
#     
#     if (length(meta2.tmp$ph[is.na(meta2.tmp$ph) == FALSE]) >= 1) {
#       grid.arrange(gg.taxa.bar + theme(legend.position="none"),
#                    g_legend(gg.taxa.bar), 
#                    gg.men.stat.dot + theme(legend.position="none"),
#                    g_legend(gg.men.stat.dot),
#                    gg.tanner.dot + theme(legend.position="none"),
#                    g_legend(gg.tanner.dot),
#                    gg.ph.line + theme(legend.position="none"),
#                    ncol=2, widths=c(4,1), heights=c(5,1,1.25,1))
#     }
#     else 
grid.arrange(gg.taxa.bar + theme(legend.position="none"),
             g_legend(gg.taxa.bar), 
             gg.men.stat.dot + theme(legend.position="none"),
             g_legend(gg.men.stat.dot),
             gg.tanner.dot + theme(legend.position="none"),
             g_legend(gg.tanner.dot),
             ncol=2, widths =c(4,1), heights=c(5,1,1.25))
  }
}
dev.off()
