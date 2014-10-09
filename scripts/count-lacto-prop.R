load("data-postproc/02-cluster-pcoa-last-run.RData")

library(reshape)
library(ggplot2)

meta$LAB <- colSums(prop.red[c(6, 28, 43:49, 72:74),]) # 6=Aerococcus, 28=Facklamia, 43:49=Lactobacillus, 72:74=Streptococcus
meta$Lacto <- colSums(prop.red[grep("Lactobacillus", rownames(prop.red)),])

## Number of girls with at least 10% LAB in vagina prior to menarche
length(unique(na.omit(meta$subject[meta$type=="girl" & meta$site=="vag" & meta$men.stat=="pre" & meta$LAB >= 0.1])))
## n=28 (118, 129 did not meet criterion; 104 had 99.8% LAB on vulva but no premenarcheal vagina swab)
## 28/31 = 90.3%

## Number of girls with at least 10% Lactobacillus in vagina prior to menarche
length(unique(na.omit(meta$subject[meta$type=="girl" & meta$site=="vag" & meta$men.stat=="pre" & meta$Lacto >= 0.1])))
## n=27 (124 missing relative to above; had Streptococcus)
## 27/31 = 87.1%%

## Number of girls with at least 50% LAB in vagina prior to menarche
length(unique(na.omit(meta$subject[meta$type=="girl" & meta$site=="vag" & meta$men.stat=="pre" & meta$LAB >= 0.5])))
## n=25 (102, 118, 124, 129, 135 did not meet criterion; N/A for 104)
## 25/31 = 80.5%

## Number of girls with at least 50% Lactobacillus in vagina prior to menarche
length(unique(na.omit(meta$subject[meta$type=="girl" & meta$site=="vag" & meta$men.stat=="pre" & meta$Lacto >= 0.5])))
## n=25 (same as above)
## 25/31 = 80.5%

meta.tmp <- meta[meta$type=="girl" & meta$site=="vag",
                 c("subject", "visit", "men.stat", "LAB", "Lacto")]

meta.tmp <- melt(meta.tmp, id.vars=c("subject", "visit", "men.stat"))

gg.girl.lab <- ggplot(meta.tmp, aes(x=visit, y=value, group=variable, color=variable, shape=men.stat))

gg.girl.lab + geom_line(lty=3, alpha=0.5) +
  facet_wrap( ~ subject, ncol=4) +
  geom_hline(yintercept=0.1, linetype=2, color="darkorange1") + 
  geom_hline(yintercept=0.5, linetype=2, color="chartreuse1") +
  geom_point(size=3) +
  scale_color_manual(values=c("gray10", "gray50"),
                     breaks=c("LAB", "Lacto"),
                     labels=c("LAB", "Lactobacillus"),
                     name="Bacteria") +
  scale_shape_manual(values=c(16, 1), 
                     breaks=c("pre", "post"), 
                     labels=c("Pre", "Post"),
                     name="Menarche\nStatus", na.value=10) +
  xlab("Visit") +
  ylab("Proportion") +
  scale_x_discrete(1:14, name="Visit No.") +
  ylim(0,1) +
  theme_cust_nogrid +
  theme(axis.text=element_text(size=6))

ggsave("girl-vag-lab-lacto.png", width=10, height=8, units="in")
 