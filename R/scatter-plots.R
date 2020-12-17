# Tissue centric scatter plot -----------------
#' Scatter plot of two columns of a data.frame after log2 transformation.
#'
#' @param data numeric data.frame
#' @param x column name to be mapped to the x-axis
#' @param y column name to be mapped to the y-axis
#' @param title plot title
#' @param xlabs x-axis title
#' @param ylabs y-axis title
#' @param min numeric, minimum of both x-axis and y-axis
#' @param max numeric, maximum of both x-axis and y-axis
#' @param a numeric, alpha scaling factor to apply to the data
#' @param xcor numeric, x coordinate for the correlation annotation to be anchored on its center
#' @param ycor numeric, y coordinate for the correlation annotation to be anchored
#' @param binx numeric, bin to use for the marginal plot of the x-axis
#' @param biny numeric, bin to use for the marginal plot of the y-axis
#' @param methodcor 'all' for both Spearman and Pearson or or 'spearman' or 'pearson'
#' @param usecor an optional character string giving a method for computing the correlation in the presence of missing values.
#'               This must be (an abbreviation of) one of the strings
#'               "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs" (default).
#' @param dens logical, default TRUE. Whether density countour curves should be added to the main figure
#' @param rug logical, default TRUE. Whether a rug should be added on the main figure
#' @param geomsize numeric, default 2. the size of the geom_point()
#' @param fig character string to chose from a set of predefined options that allows different outputs
#'           'all' : complete figure with annotation and marginal plots
#'           'grob': list of the grobs use by gridExtra::arrangeGrob to create the plot
#'           'list': list of the different plots
#'           'p'   : the main figure with the legend
#'           'p1'  : the main figure without the legend
#'           'p2'  : the marginal figure on the x-axis
#'           'p3'  : the marginal figure on the y-axis
#'           'legend': the legend of the figure only
#' @param scaleCol colour palette to use for the plot
#' @param verbose logical; default FALSE
#' @param smooth logical; default TRUE whether a smooth curve should be plotted
#' @param smoothMethod string, method to use for the smoothing default: "gam" (from mgcv package)
#' @param smoothFormula formula to use with the smooth method default: as.formula("y ~ s(x, bs = 'cs')")
#' @param xmin numeric, minimum of the x-axis
#' @param ymin numeric, minimum of the y-axis
#' @param xmax numeric, maximum of the x-axis
#' @param ymax numeric, maximum of the y-axis
#' @param abline logical; default: TRUE; whether the line y=x should be added as a reference
#' @param publi logical, default: TRUE; whether to apply the theme_bw which is better suited for publication
#' @param col colour palette
#' @param labCol title for the colour scale
#' @param family police family to use for the labels and titles
#' @param policeSize base size for the plot
#' @param sizeLegend size of the legend
#' @param sizeCor size of the correlation annotation; default: 3
#' @param dealWithSpace logical; default FALSE.
#'                      Whether the space in the names of the columns to be plotted should be removed.
#' @param p2alpha alpha scale for the marginal plot on the x-axis
#' @param p3alpha alpha scale for the marginal plot on the y-axis
#'
#' @return a plot (see parameter fig for better description)
#'
#' @importFrom mgcv gam
#' @export
#'
scatplot.log2<-function(data,x,y,title="",xlabs,ylabs,min,max,a=0.65,
                        xcor=3.3,ycor=16.999,binx=0.2,biny=0.2,
                        methodcor='all',usecor="pairwise.complete.obs",
                        dens=TRUE,rug=TRUE,geomsize=2,
                        fig='all',scaleCol,verbose=FALSE,
                        smooth=TRUE, smoothMethod="gam", smoothFormula=as.formula("y ~ s(x, bs = 'cs')"),
                        xmin,ymin,xmax,ymax,abline=TRUE,publi=TRUE,col=NULL,labCol,
                        family='',policeSize=11, sizeLegend=8, sizeCor=3, dealWithSpace=FALSE,
                        p2alpha=1,p3alpha=1){

  if(dealWithSpace) {
    colnames(data)<-gsub(" ","",colnames(data))
    x<-gsub(' ','',x)
    y<-gsub(' ','',y)
    colnames(data)<-gsub('\\(','',colnames(data))
    x<-gsub('\\(','',x)
    y<-gsub('\\(','',y)
    colnames(data)<-gsub('\\)','',colnames(data))
    x<-gsub('\\)','',x)
    y<-gsub('\\)','',y)
  }

  ##creation of the main plot
  x2=paste(paste("log2(",substitute(x)),")")
  y2=paste(paste("log2(",substitute(y)),")")

  if(missing(xlabs)) xlabs=paste(paste("log2(",substitute(x)),")")
  if(missing(ylabs)) ylabs=paste(paste("log2(",substitute(y)),")")

  BW<-TRUE
  if(!is.null(col)){
    col<-substitute(col)
    BW<-FALSE
    if(missing(labCol))
      labCol<-col
  }

  if(missing(xmax)){
    if(missing(max)){
      xmax<-base::max(log2(data[,x]+1))+1
    }else{
      xmax<-max
    }
  }
  if(missing(ymax)){
    if (missing(max)) {
      ymax <- base::max(log2(data[,y]+1))+1
    }else{
      ymax<-max
    }
  }
  if(missing(xmin)){
    if (missing(min)){
      xmin<-base::min(log2(data[,x]+1))-1
    }else{
      xmin <- min
    }
  }

  if(missing(ymin)){
    if (missing(min)){
      ymin<-base::min(log2(data[,y]+1))-1
    }else{
      ymin<-min
    }
  }

  row_sub <- apply(data[,c(x,y)], 1, function(row) all(row !=0 ))
  data2   <- data[row_sub,]
  if(verbose)   print(paste(nrow(data)-nrow(data2), ' points were not plotted since the value was null in one case'))

  if (length(col)>0){
    p <- ggplot(data2, aes_string(x=x2, y=y2, colour=col))
    if(!missing(scaleCol)){
      p <- p + scale_colour_manual(values=scaleCol)
    }
    p <- p + geom_point(alpha=a)
  }else{
    p <- ggplot(data2, aes_string(x=x2, y=y2))
    p <- p+geom_point(alpha=a)
  }

  p <- p + labs(title=title, x=xlabs, y=ylabs)
  p <- p + coord_cartesian(xlim=c(xmin,xmax),ylim=c(ymin,ymax))

  if (abline) p <- p + geom_abline(colour = "black")

  #for the r coefficient (on the graph):
  switch(methodcor,
         'all'={p <- p + annotate('text',
                                  family=family,
                                  size=sizeCor,
                                  label=paste("\u03C1 (Spearman) =",
                                              signif(cor(log2(data[x]+1),
                                                         log2(data[y]+1),
                                                         use=usecor,
                                                         method='spearman')[1,1],
                                                     digits=3),
                                              " ; ",
                                              "r (Pearson) =",
                                              signif(cor(log2(data[x]+1),
                                                         log2(data[y]+1),
                                                         use=usecor,
                                                         method='pearson')[1,1],
                                                     digits=3)),
                                  x=xcor,y=ycor)},
         {p <- p + annotate('text',
                            family=family,
                            size=sizeCor,
                            label=paste0("r (",
                                         simpleCap(methodcor),
                                         ") = ",
                                         signif(cor(log2(data[x]+1),
                                                    log2(data[y]+1),
                                                    use=usecor,
                                                    method=methodcor)[1,1],
                                                digits=3)),
                            x=xcor,y=ycor)
         }
  )


  if(smooth)  p <- p + geom_smooth(aes_string(x=x2,y=y2),
                                   data=data2,
                                   method=smoothMethod,
                                   formula=smoothFormula
  )
  if(dens){#contour (circles are circonvening  area of same points density)
    p <- p + geom_density2d(data=data2,na.rm=TRUE,aes_string(x=x2,y=y2,alpha='..level..', colour=NULL))
  }

  if (publi) {
    p <- p + theme_bw(base_family=family, base_size=policeSize)
    p <- p + theme(panel.border=element_rect(linetype=NULL,colour=NA))
    p <- p + theme(plot.title=element_text(hjust=0.5))
  }

  #the legend is formatted and then it is extracted from the main plot and would be added at the end


  if(!BW){
    p <- p + guides(ncol=2,
                    colour = guide_legend(title=labCol,
                                          keywidth=0.5,keyheight=0.5,title.position = "right",
                                          title.theme = element_text(face='bold',size=sizeLegend,angle = 90,family = family ),
                                          label.theme = element_text(size=sizeLegend,angle=0,family = family)),
                    alpha = guide_legend(keywidth=0.5,keyheight=0.5,title.position = "right",
                                         title.theme = element_text(face='bold',size=sizeLegend,angle = 90,family=family),
                                         label.theme = element_text(size=sizeLegend,angle=0,family = family))
    )
  }else{
    p<- p +guides(ncol=2,
                  alpha=guide_legend(keywidth=0.5,keyheight=0.5,title.position = "top",
                                     title.theme = element_text(face='bold',size=sizeLegend,angle = 0,family=family),
                                     label.theme = element_text(size=sizeLegend,angle=0,family=family)))
  }

  legend <- g_legend(p)

  if (rug){#allows some sort of visualization of the most dense parties of one axis
    if(BW){
      p <- p + geom_rug(col=rgb(0,0,0.5,alpha=0.015))
    }else{
      p <- p + geom_rug(alpha=0.25)
    }
  }
  p1 <- p + theme(legend.position='none')#legend stripped from this plot

  ##creation of the histogram for the x axis
  #note: other possible way: instead of count => density
  p2 <- ggplot(data2, aes_string(x=x2, fill=col,alpha=p2alpha))
  if (!missing(scaleCol))
    p2 <- p2 + scale_fill_manual(values=scaleCol)
  if(publi)   p2<- p2 +theme_bw(base_family=family, base_size=policeSize) +theme(panel.border=element_rect(linetype=NULL,colour=NA))

  p2 <- p2 + geom_histogram(aes(y=..count..),binwidth=binx)
  p2 <- p2 + coord_cartesian(xlim=c(xmin,xmax))
  p2 <- p2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
  p2 <- p2 + scale_y_reverse()
  p2 <- p2 + theme(axis.text.x=element_blank())
  p2 <- p2 + theme(legend.title = element_blank(),legend.position='none')

  ##creation of the histogram for the y axis
  #note other possible way: instead of count => density
  p3 <- ggplot(data2, aes_string(x=y2, fill=col,alpha=p3alpha))
  if (!missing(scaleCol))
    p3 <- p3+  scale_fill_manual(values=scaleCol)
  if(publi)   p3<-p3+theme_bw(base_family=family, base_size=policeSize) +theme(panel.border=element_rect(linetype=NULL,colour=NA))

  p3 <- p3 + geom_histogram(aes(y=..count..),binwidth=biny)
  p3 <- p3 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
  p3 <- p3 + scale_y_reverse()
  p3 <- p3 + theme(axis.text.x  = element_text(angle=90, vjust=0,hjust=1))
  p3 <- p3 + coord_flip(xlim=c(ymin,ymax) )
  p3 <- p3 + theme(axis.text.y=element_blank())
  p3 <- p3 + theme(legend.position="none")

  #creation of the object before drawing them
  gp1<-ggplot_gtable(ggplot_build(p1))
  gp2<-ggplot_gtable(ggplot_build(p2))
  gp3<-ggplot_gtable(ggplot_build(p3))

  #to synchronize the x axis of the main plot with the histogram of the x variable
  gp2$widths <- gp1$widths
  #to synchronize the x axis of the main plot with the histogram of the y variable
  gp3$heights <-gp1$heights

  switch(fig,
         'all' =   {return(grid.arrange(arrangeGrob(gp3, gp1,legend,gp2,  widths=c(1,5), heights=c(5,1))))},
         'grob'=   {return(arrangeGrob(gp3, gp1,legend,gp2,  widths=c(1,5), heights=c(5,1)))},
         'list'=   {return(list(gp3, gp1,legend,gp2))},#changed to list of something
         'p'   =   {return(print(p))},
         'p1'  =   {return(print(p1))},
         'p2'  =   {return(print(p2))},
         'p3'  =   {return(print(p3))},
         'legend'= {return(print(legend))},
         {return(p)}
  )
}

#' Scatter plot of two columns of a data.frame after possible log2 transformation.
#'
#' @param data numeric data.frame
#' @param x column name to be mapped to the x-axis
#' @param y column name to be mapped to the y-axis
#' @param title plot title
#' @param xlabs x-axis title
#' @param ylabs y-axis title
#' @param min numeric, minimum of both x-axis and y-axis
#' @param max numeric, maximum of both x-axis and y-axis
#' @param a numeric, alpha scaling factor to apply to the data
#' @param log2 logical, default TRUE. Whether the data should be log2 transformed before plotting
#' @param xcor numeric, x coordinate for the correlation annotation to be anchored on its center
#' @param ycor numeric, y coordinate for the correlation annotation to be anchored
#' @param binx numeric, bin to use for the marginal plot of the x-axis
#' @param biny numeric, bin to use for the marginal plot of the y-axis
#' @param methodcor 'all' for both Spearman and Pearson or or 'spearman' or 'pearson'
#' @param usecor an optional character string giving a method for computing the correlation in the presence of missing values.
#'               This must be (an abbreviation of) one of the strings
#'               "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs" (default).
#' @param dens logical, default TRUE. Whether density countour curves should be added to the main figure
#' @param rug logical, default TRUE. Whether a rug should be added on the main figure
#' @param geomsize numeric, default 2. the size of the geom_point()
#' @param fig character string to chose from a set of predefined options that allows different outputs
#'           'all' : complete figure with annotation and marginal plots
#'           'grob': list of the grobs use by gridExtra::arrangeGrob to create the plot
#'           'list': list of the different plots
#'           'p'   : the main figure with the legend
#'           'p1'  : the main figure without the legend
#'           'p2'  : the marginal figure on the x-axis
#'           'p3'  : the marginal figure on the y-axis
#'           'legend': the legend of the figure only
#' @param scaleCol colour palette to use for the plot
#' @param smooth logical; default TRUE whether a smooth curve should be plotted
#' @param smoothMethod string, method to use for the smoothing default: "gam" (from mgcv package)
#' @param smoothFormula formula to use with the smooth method default: as.formula("y ~ s(x, bs = 'cs')")
#' @param xmin numeric, minimum of the x-axis
#' @param ymin numeric, minimum of the y-axis
#' @param xmax numeric, maximum of the x-axis
#' @param ymax numeric, maximum of the y-axis
#' @param abline logical; default: TRUE; whether the line y=x should be added as a reference
#' @param publi logical, default: TRUE; whether to apply the theme_bw which is better suited for publication
#' @param col colour palette
#' @param labCol title for the colour scale
#' @param verbose logical; default FALSE
#' @return a plot (see parameter fig for better description)
#' @export
#'
scatplot<-function(data,x,y,title,xlabs,ylabs,min,max,a=0.65,log2=TRUE,
                   xcor=3.3,ycor=16.999,binx=0.2,biny=0.2,
                   methodcor='all',usecor="pairwise.complete.obs",
                   dens=TRUE,rug=TRUE,geomsize=2,fig='all',scaleCol,
                   smooth=TRUE,smoothMethod="gam", smoothFormula=as.formula("y ~ s(x, bs = 'cs')"),
                   xmin,ymin,xmax,ymax,abline=TRUE,publi=TRUE,col=NULL,labCol,
                   verbose=FALSE){

  ##creation of the main plot
  if(log2){
    x2=paste(paste("log2(",substitute(x)),")")
    y2=paste(paste("log2(",substitute(y)),")")
  }else{
    x2=substitute(x)
    y2=substitute(y)
  }


  BW<-TRUE
  if(!is.null(col)){
    col<-substitute(col)
    BW<-FALSE
    if(missing(labCol))
      labCol<-col
  }

  if(missing(xmax)){
    if(missing(max)){
      xmax <- base::max(x,y)
    }else{
      xmax<-max
    }
  }
  if(missing(ymax)){
    if (missing(max)) {
      ymax <- base::max(y,x)
    }else{
      ymax<-max
    }
  }
  if(missing(xmin)){
    if (missing(min)){
      xmin <- base::min(x,y)
    }else{
      xmin <- min
    }
  }

  if(missing(ymin)){
    if (missing(min)){
      ymin<- base::min(y,x)
    }else{
      ymin<-min
    }
  }
  if(log2){
    row_sub <- apply(data[,c(x,y)], 1, function(row) all(row !=0 ))
    data2<-data[row_sub,]
    if(verbose) print(paste(nrow(data)-nrow(data2), ' points were not plotted since the value was null in one case'))
  }else{
    data2<-data
  }
  if (length(col)>0){
    p<-ggplot(data2, aes_string(x=x2, y=y2, colour=col))
    if(!missing(scaleCol)){
      p <- p + scale_colour_manual(values=scaleCol)
    }
    p<- p+geom_point(alpha=a)
  }else{
    p <- ggplot(data2, aes_string(x=x2, y=y2))
    p <- p+geom_point(alpha=a)
  }
  if (publi) {
    p<- p + theme_bw() + theme(panel.border=element_rect(linetype=NULL,colour=NA))
  }

  p <- p + labs(title=title, x=xlabs, y=ylabs)
  p <- p + coord_cartesian(xlim=c(xmin,xmax),ylim=c(ymin,ymax))

  if (abline) p <- p + geom_abline(colour = "black")

  #for the r coefficient (on the graph):
  if(log2){
    switch(methodcor,
           'all'={p <-p + annotate('text',
                                   label=paste("\u03C1 (Spearman)",
                                               signif(cor(log2(data[x]+1),
                                                          log2(data[y]+1),
                                                          use=usecor,
                                                          method='spearman')[1,1],
                                                      digits=3),
                                               " ; ",
                                               "r (Pearson) =",
                                               signif(cor(log2(data[x]+1),
                                                          log2(data[y]+1),
                                                          use=usecor,
                                                          method='pearson')[1,1],
                                                      digits=3)),
                                   x=xcor,
                                   y=ycor)
           },
           {p <- p + annotate('text',
                              label=paste0("r (",
                                           simpleCap(methodcor),
                                           ") = ",
                                           signif(cor(log2(data[x]+1),
                                                      log2(data[y]+1),
                                                      use=usecor,
                                                      method=methodcor)[1,1],
                                                  digits=3)),
                              x=xcor,
                              y=ycor)
           }
    )
  }else{
    switch(methodcor,
           'all'={p <-p + annotate('text',
                                   label=paste("\u03C1 (Spearman)",
                                               signif(cor(data[x],
                                                          data[y],
                                                          use=usecor,
                                                                  method='spearman')[1,1],
                                                      digits=3),
                                               " ; ",
                                               "r (Pearson) =",
                                               signif(cor(data[x],
                                                          data[y],
                                                          use=usecor,
                                                          method='pearson')[1,1],
                                                      digits=3)),
                                   x=xcor,
                                   y=ycor)
           },
           {p <-p + annotate('text',
                             label=paste0("r (",
                                          simpleCap(methodcor),
                                          ") = ",
                                          signif(cor(data[x],
                                                     data[y],
                                                     use=usecor,
                                                     method=methodcor)[1,1],
                                                 digits=3)),
                             x=xcor,
                             y=ycor)
           }
    )
  }


  if(smooth)  p <- p + geom_smooth(aes_string(x=x2,y=y2),
                                   data=data2,#data[data[,x] > 0 &  data[,y]>0,],
                                   method=smoothMethod,
                                   formula=smoothFormula
  )
  if(dens){#contour (circles are circonvening  area of same points density)
    p <- p + geom_density2d(data=data2,na.rm=TRUE,aes_string(x=x2,y=y2,alpha='..level..', colour=NULL))
  }
  #the legend is formatted and then it is extracted from the main plot and would be added at the end

  if(!BW){
    p <- p + guides(ncol=2,
                    colour = guide_legend(title=labCol,
                                          keywidth=0.5,keyheight=0.5,title.position = "right",
                                          title.theme = element_text(face='bold',size=8,angle = 90),
                                          label.theme = element_text(size=8,angle=0)),
                    alpha = guide_legend(keywidth=0.5,keyheight=0.5,title.position = "right",
                                         title.theme = element_text(face='bold',size=8,angle = 90),
                                         label.theme = element_text(size=8,angle=0))
    )
  }else{
    p<- p +guides(ncol=2,
                  alpha=guide_legend(keywidth=0.5,keyheight=0.5,title.position = "right",
                                     title.theme = element_text(face='bold',size=8,angle = 90),
                                     label.theme = element_text(size=8,angle=0)))
  }

  legend <- g_legend(p)

  if (rug){#allows some sort of visualization of the most dense parties of one axis
    if(BW){
      p <- p + geom_rug(col=rgb(0,0,0.5,alpha=0.015))
    }else{
      p <- p + geom_rug(alpha=0.25)
    }
  }
  p1 <- p + theme(legend.position='none')#legend stripped from this plot

  ##creation of the histogram for the x axis
  #note: other possible way: instead of count => density
  p2 <- ggplot(data2, aes_string(x=x2, fill=col))
  if (!missing(scaleCol))
    p2 <- p2 + scale_fill_manual(values=scaleCol)
  if(publi)   p2<- p2 +theme_bw() +theme(panel.border=element_rect(linetype=NULL,colour=NA))

  p2 <- p2 + geom_histogram(aes(y=..count..),binwidth=binx)
  p2 <- p2 + coord_cartesian(xlim=c(xmin,xmax))
  p2 <- p2 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
  p2 <- p2 + scale_y_reverse()
  p2 <- p2 + theme(axis.text.x=element_blank())
  p2 <- p2 + theme(legend.title = element_blank(),legend.position='none')

  ##creation of the histogram for the y axis
  #note other possible way: instead of count => density
  p3 <- ggplot(data2, aes_string(x=y2, fill=col))
  if (!missing(scaleCol))
    p3 <- p3+  scale_fill_manual(values=scaleCol)
  if(publi)   p3<-p3+theme_bw() +theme(panel.border=element_rect(linetype=NULL,colour=NA))

  p3 <- p3 + geom_histogram(aes(y=..count..),binwidth=biny)
  p3 <- p3 + theme(axis.title.y=element_blank(), axis.title.x=element_blank())
  p3 <- p3 + scale_y_reverse()
  p3 <- p3 + theme(axis.text.x  = element_text(angle=45, vjust=0,hjust=1))
  p3 <- p3 + coord_flip(xlim=c(ymin,ymax) )
  p3 <- p3 + theme(axis.text.y=element_blank())
  p3 <- p3 + theme(legend.position="none")

  #creation of the object before drawing them
  gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
  gp2<-suppressWarnings(ggplot_gtable(ggplot_build(p2)))
  gp3<-suppressWarnings(ggplot_gtable(ggplot_build(p3)))

  #to synchronize the x axis of the main plot with the histogram of the x variable
  gp2$widths <- gp1$widths
  #to synchronize the x axis of the main plot with the histogram of the y variable
  gp3$heights <-gp1$heights

  switch(fig,
         'all' =   {return(suppressWarnings(grid.arrange(arrangeGrob(gp3, gp1,legend,gp2,  widths=c(1,5), heights=c(5,1)))))},
         'grob'=   {return(suppressWarnings(arrangeGrob(gp3, gp1,legend,gp2,  widths=c(1,5), heights=c(5,1))))},
         'list'=   {return(list(gp3, gp1,legend,gp2))},#changed to list of something
         'p'   =   {return(print(p))},
         'p1'  =   {return(print(p1))},
         'p2'  =   {return(print(p2))},
         'p3'  =   {return(print(p3))},
         'legend'= {return(print(legend))},
         {return(p)}
  )
}




# Gene centric scatter plot -----------------------------

#' Scatter plot for the expression of a gene for two datasets across different tissues
#'
#' @param ID character string, gene identifier
#' @param DF1 numeric data.frame or equivalent
#' @param DF2 numeric data.frame or equivalent
#' @param palette named character vector with the colour hex code to use
#' @param nameIt logical, default TRUE.
#'               Whether the gene name (from mapID)
#'               should be incorporated to the plot title
#' @param mapID named character vector. Map the gene identifiers to their name.
#'              default: \code{\link[barzinePhdData:gene.mapID]{gene.mapID}}
#' @param logIt logical, default TRUE. Whether the data should be log2 transform before plotting
#' @param label.DF1 character string for DF1 labelling. Used as x-axis title
#' @param label.DF2 character string for DF2 labelling. Used as y-axis title
#' @param xcoord vector of two numeric. Gives the min and max of the x-axis
#' @param ycoord vector of two numeric. Gives the min and max of the y-axis
#' @param colNamesDF1 default: c('Tissue','Transcriptomics')
#' @param colNamesDF2 default: c('Tissue','Proteomics')
#' @param base_size numeric, default: 11
#' @param base_family default: Linux Libertine
#' @param label_size numeric, default 3.4
#' @param label_family default: same as base_family
#' @param centerTitle logical, default TRUE
#' @param useCor default: 'everything'
#'               an optional character string for 'cor'
#'               which gives a method for computing covariances
#'               in the presence of missing values.
#'               This must be (an abbreviation of) one of the strings
#'               "everything", "all.obs", "complete.obs",
#'               "na.or.complete", or "pairwise.complete.obs".
#'
#' @return a plot
#' @export
#'
geneScatters<-function(ID,DF1,DF2,palette,
                       nameIt=TRUE,mapID=barzinePhdData::gene.mapID,
                       logIt=TRUE,
                       label.DF1="Uhl\u00E9n et al. (mRNA)\nlog2(FPKM+1)",
                       label.DF2='Pandey lab (protein)\nlog2(PPKM+1)',
                       xcoord,ycoord,
                       colNamesDF1=c('Tissue','Transcriptomics'),
                       colNamesDF2=c(colNamesDF1[1],'Proteomics'),
                       base_size=11,base_family='Linux Libertine',
                       label_size=3.4,label_family=base_family,
                       centerTitle=TRUE,
                       useCor='everything'
){

  rho =expression(rho)
  t1<-reshape2::melt(DF1[ID,])
  t2<-reshape2::melt(DF2[ID,])
  colnames(t1)<-colNamesDF1
  colnames(t2)<-colNamesDF2

  newDF<-merge(t1,t2,by=colNamesDF1[1])
  if(logIt){
    p <- ggplot(newDF, aes_string(x=paste0('log2(',colNamesDF1[2],'+1)'),
                                  y=paste0('log2(',colNamesDF2[2],'+1)'),
                                  color=colNamesDF1[1],label=colNamesDF1[1]))
  }else{
    p <- ggplot(newDF, aes_string(x=colNamesDF1[2],y=colNamesDF2[2],
                                  color=colNamesDF1[1],label=colNamesDF1[1]))
  }
  p <- p + geom_point()
  p <- p + scale_colour_manual(values=palette)
  p <- p + geom_text_repel(size=label_size,family=label_family)
  p <- p + theme_bw(base_size=base_size, base_family = base_family)
  p <- p + theme(legend.position = 'none')
  if(nameIt){
    Title=paste0(mapID[ID],' (',ID,')\n')
  }else{
    Title=paste0(ID,'\n')
  }
  if(logIt){
    Title=paste0(Title,
                 'Spearman \u03C1 = ',
                 signif(cor(x=log2(as.numeric(DF1[ID,])+1),
                            y=log2(as.numeric(DF2[ID,]+1)),
                            method='spearman',
                            use=useCor),digits = 2),' \u2022 ',
                 'Pearson r = ',
                 signif(cor(x=log2(as.numeric(DF1[ID,])+1),
                            y=log2(as.numeric(DF2[ID,])+1),
                            use=useCor),digits = 2))
  }else{
    Title=paste0(Title,
                 '\u03C1 (Spearman) = ',
                 signif(cor(x=as.numeric(DF1[ID,]),
                            y=as.numeric(DF2[ID,]),
                            method='spearman',
                            use=useCor),digits = 2),' \u2022 ',
                 'r (Pearson) = ',
                 signif(cor(x=as.numeric(DF1[ID,]),
                            y=as.numeric(DF2[ID,]),
                            use=useCor),digits = 2))
  }
  p <- p + labs(x=label.DF1,y=label.DF2,title=Title)
  if(!missing(xcoord)) p<-p +coord_cartesian(xlim=xcoord)
  if(!missing(ycoord)) p<-p +coord_cartesian(ylim=xcoord)

  if(centerTitle) p <- p +theme(plot.title=element_text(hjust=0.5))

  return(p)
}


utils::globalVariables("count")
