library(RColorBrewer)

df <- read.csv("data/79taxa.csv", header=TRUE)

## There are 17 unique orders and 7 unique phyla

YlGnBu <- brewer.pal(9, "YlGnBu")
YlOrRd <- brewer.pal(9, "YlOrRd")
RdPu <- brewer.pal(9, "RdPu")
YlGn <- brewer.pal(9, "YlGn")
PuOr <- brewer.pal(11, "PuOr")
Greys <- brewer.pal(9, "Greys")
BrBG <- brewer.pal(11, "BrBG")

df$order.col <- rep("NA", nrow(df))

lacto.col <- colorRampPalette(c(YlGnBu[2:4], "cyan3", "dodgerblue2", YlGnBu[7:9]))(length(df$order.col[grep("Lactobacillales", df$Order)]))
df$order.col[grep("Lactobacillales", df$Order)] <- lacto.col[c(1,2,4,3,5,6,7,9,8,10,11,12)]
# pie(rep(1,12), col=df$order.col[grep("Lactobacillales", df$Order)], labels=df$taxon[grep("Lactobacillales", df$Order)])

df$order.col[grep("Clostridiales", df$Order)] <- colorRampPalette(c(YlOrRd, BrBG[1:4]))(length(df$order.col[grep("Clostridiales", df$Order)]))
# pie(rep(1,23), col=df$order.col[grep("Clostridiales", df$Order)], labels=df$taxon[grep("Clostridiales", df$Order)])

df$order.col[grep("Actinomycetales", df$Order)] <- colorRampPalette(c("lightgoldenrod1", "olivedrab1", YlGn[c(7,9)]))(length(df$order.col[grep("Actinomycetales", df$Order)]))
# pie(rep(1,9), col=df$order.col[grep("Actinomycetales", df$Order)], labels=df$taxon[grep("Actinomycetales", df$Order)])

df$order.col[grep("Bifidobacteriales", df$Order)] <- RdPu[5:8]
# pie(rep(1,4), col=df$order.col[grep("Bifidobacteriales", df$Order)], labels=df$taxon[grep("Bifidobacteriales", df$Order)])

df$order.col[grep("Coriobacteriales", df$Order)] <- PuOr[c(10,11,8,9)]
# pie(rep(1,4), col=df$order.col[grep("Coriobacteriales", df$Order)], labels=df$taxon[grep("Coriobacteriales", df$Order)])

df$order.col[grep("Bacteroidales", df$Order)] <- BrBG[7:11]
# pie(rep(1,5), col=df$order.col[grep("Bacteroidales", df$Order)], labels=df$taxon[grep("Bacteroidales", df$Order)])

df$order.col[!(df$Order %in% c("Lactobacillales",
                               "Clostridiales",
                               "Actinomycetales",
                               "Bifidobacteriales",
                               "Coriobacterialies",
                               "Bacteroidales"))] <- colorRampPalette(Greys[2:8])(length(df$order.col[!(df$Order %in% c("Lactobacillales",
                                                                                                                       "Clostridiales",
                                                                                                                       "Actinomycetales",
                                                                                                                       "Bifidobacteriales",
                                                                                                                       "Coriobacterialies",
                                                                                                                       "Bacteroidales"))]))

## Write new taxa color palette
col.taxa <- df$order.col
names(col.taxa) <- df$taxon

## Color only LAB spp. -- this was for a presentation I made
# df$col.LAB <- rep("#F2F2F2", nrow(df))
# df$col.LAB[grep("Lactobacillales", df$Order)] <- lacto.col[c(1,2,4,3,5,6,7,9,8,10,11,12)]
# col.LAB <- df$col.LAB
# names(col.LAB) <- df$taxon

## Color only Gardnerella -- again, just for a presentation
# df$col.GV <- rep("#F2F2F2", nrow(df))
# df$col.GV[grep("Gardnerella", df$Genus)] <- RdPu[7:8]
# col.GV <- df$col.GV
# names(col.GV) <- df$taxon