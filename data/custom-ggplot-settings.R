## Custom graphing options

## Custom ggplot theme to call when you make plots

library(ggplot2)

theme_cust <- theme_bw() +
  theme(plot.title=element_text(face="bold", size=14),
        axis.title=element_text(size=10),
        axis.text=element_text(size=8),
        legend.title=element_text(face="bold", size=10),
        legend.text=element_text(size=8),
        legend.background=element_blank(),
        strip.text=element_text(face="bold", size=8), 
        strip.background=element_blank(),
        panel.background=element_rect(fill="transparent", color = NA),
        plot.background=element_rect(fill="transparent", color = NA))

## omit minor gridlines
theme_cust_nominor <- theme_cust +
  theme(panel.grid.minor=element_blank())

## omit major and minor gridlines
theme_cust_nogrid <- theme_cust +
  theme(panel.grid.minor=element_blank(),
        panel.grid.major=element_blank())

# Multiple plot function with layout option from http://www.cookbook-r.com/Graphs/Multiple_graphs_on_one_page_(ggplot2)/
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

## Handy function for making semi-transparent colors
makeTransparent = function(..., alpha=0.5) {
  if(alpha<0 | alpha>1) stop("alpha must be between 0 and 1")
  
  alpha = floor(255*alpha)  
  newColor = col2rgb(col=unlist(list(...)), alpha=FALSE)
  
  .makeTransparent = function(col, alpha) {
    rgb(red=col[1], green=col[2], blue=col[3], alpha=alpha, maxColorValue=255)
  }
  
  newColor = apply(newColor, 2, .makeTransparent, alpha=alpha)
  
  return(newColor)
  
} # thanks to http://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color