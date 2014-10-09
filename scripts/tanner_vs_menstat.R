## Determine proportions of pre/post in each Tanner stage

## Tanner breast
tb.ms <- matrix(nrow=5, ncol=2)
for (i in c(1,2,3,4,5)){
  tb.ms[i,1] <- length(na.exclude(meta.vag$tan.br.dr[meta.vag$type=="girl" & meta.vag$men.stat=="pre" & meta.vag$tan.br.dr==i]))
  tb.ms[i,2] <- length(na.exclude(meta.vag$tan.br.dr[meta.vag$type=="girl" & meta.vag$men.stat=="post" & meta.vag$tan.br.dr==i]))
}
rownames(tb.ms) <- c("TB-1","TB-2","TB-3","TB-4","TB-5")
colnames(tb.ms) <- c("pre","post")

## Tanner genital
tg.ms <- matrix(nrow=5, ncol=2)
for (i in c(1,2,3,4,5)){
  tg.ms[i,1] <- length(na.exclude(meta.vag$tan.br.dr[meta.vag$type=="girl" & meta.vag$men.stat=="pre" & meta.vag$tan.gen.dr==i]))
  tg.ms[i,2] <- length(na.exclude(meta.vag$tan.br.dr[meta.vag$type=="girl" & meta.vag$men.stat=="post" & meta.vag$tan.gen.dr==i]))
}
rownames(tg.ms) <- c("TG-1","TG-2","TG-3","TG-4","TG-5")
colnames(tg.ms) <- c("pre","post")

## Combine
tbg.ms <- rbind(tb.ms, tg.ms)
tbg.ms.prop <- prop.table(tbg.ms, 1)
site <- c(rep("Breast",5), rep("Genital",5))
tanner <- c(rep(c(1,2,3,4,5), 2))
tbg.ms.site.prop <- data.frame(cbind(tbg.ms.prop, site, tanner))
tbg.ms.site.prop$pre <- as.numeric(unfactor(tbg.ms.site.prop$pre))
tbg.ms.site.prop$post <- as.numeric(unfactor(tbg.ms.site.prop$post))

library(reshape)
tbg.ms.site.prop.lg <- melt(tbg.ms.site.prop, id.vars=c("site", "tanner"))

library(ggplot2)
gg.tbg.ms <- ggplot(subset(tbg.ms.site.prop.lg, tanner %in% c(2,3,4,5)), 
                           aes(x=tanner, y=value, fill=variable))

gg.tbg.ms + geom_bar(stat="identity") +
  facet_wrap(~ site, nrow=2) +
  scale_fill_manual(values=col.men.stat, 
                    name="Menarche\nstatus") +
  xlab("Tanner stage") +
  ylab("Proportion of samples") +
  theme_cust +
  theme(legend.text=element_text(size=6),
        legend.title=element_text(size=8),
        legend.position="bottom")

ggsave("tanner-vs-menstat.pdf", width=4, height=6, units="in")
