# Tissue oriented ----------------

#' Present the number of genes that are expressed in each tissue
#'
#' @param DF numeric data.frame that contains expresion data
#' @param reorder boolean. Default: TRUE. Whether the tissues
#' should be sorted in decreasing order of expressed genes
#' @param xlab character string. Label for the x-axis
#' @param ylab character string. Label for the y-axis
#' @param title character string. Title on the figure
#' @param Y a length-2 numeric vector. Limites of the y-axis
#' @param thm function theme( ) with the list of additional elements to apply
#' @param publish boolean. Default: TRUE.
#'
#' @return a figure
#' @export
#'
nbGenesPerCond<-function(DF,reorder=TRUE, xlab,ylab,title,Y,thm,
                         publish=TRUE){
  colNames<-colnames(DF) #Genes.ID are going to be added as a column
  names(colNames)<-colnames(DF) #so tempo would be a named list

  DF$Genes.ID<-rownames(DF)

  tempo<-lapply(colNames,function (x){
    VEC<-DF[,x]
    names(VEC)<-rownames(DF)
    VEC<-VEC[VEC>0]
  }
  )

  DF<-reshape2::melt(DF)
  DF<-DF[DF$value>0,]
  DF$unique<-"NotUnique"

  DF<-Reduce(rbind,lapply(colNames,function (x) {
    DIFF<-setdiff(colNames,x) #all the other conditions that will be tested against

    #next line gives the names that are unique to the current condition x
    VEC<-setdiff(names(tempo[[x]]),
                 unique(Reduce(c,lapply(DIFF,function(x) names(tempo[[x]])))) )
    if( dim(DF[DF$Genes.ID %in% VEC &DF$variable==x,])[1]>0){
      DF[DF$Genes.ID %in% VEC &DF$variable==x,]$unique<-'Unique'
    }
    DF[DF$variable==x,]
  }))

  DF$unique<-as.factor(DF$unique)
  DF$unique<-relevel(DF$unique,"Unique")

  if(reorder){
    v.order<-sapply(colNames,function(x) length(tempo[[x]]))
    #so the conditions are sorted by number of expressed genes.
    colNames<-names(v.order[order(v.order,decreasing=TRUE)])
    DF$variable<-factor(DF$variable,levels=colNames)#by
  }

  p<-ggplot(DF, aes(x=variable))+geom_bar(aes(y=..count..,fill=unique))
  # if(annot==TRUE)
  #     p<-p+stat_bin(geom="text", aes(label=..count..), vjust=-0.5, position="identity",size=4)
  if(!missing(xlab)) p <- p + xlab(xlab)
  if(!missing(ylab)) p <- p + ylab(ylab)
  if(!missing(title)) p <- p + ggtitle(title)
  if(publish) p <- p + theme_bw()
  if(!missing(Y)) p <- p + Y
  p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
  if(!missing(thm)) p <-p + thm
  print(p)

  invisible(DF)
}


#' Plot the distribution of genes across the columns.
#' Genes found uniquely in one condition are highlighted in another colour.
#'
#' @param DF numerical data.frame. Rownames have to be meaningful.
#' @param reorder boolean; default: TRUE.
#'                Whether the columns should be reorder by the number of expressed genes instead of their names.
#' @param xlab character string; used as x-axis label.
#' @param ylab character string; used as y-axis label.
#' @param title character string; title of the plot.
#' @param legend.title character string; what to use for the colour legend
#' @param publish boolean; default: TRUE. Applies \code{\link[ggplot2]{theme_bw}} to the plot.
#' @param threshold numeric; default: 0. Minimal level of expression to be considered as an expressed gene
#' @param thm other theme object (\code{\link[ggplot2]{theme}}) to pass on to the plot
#' @param Y numeric vector of two elements. Specify the limites of the y-axis.
#'          'c(minimal limit value, maximal limit value)'
#' @param remove.legend boolean; default: TRUE. Whether the legend should be removed
#'                      (do not impact the plot)
#' @param legendPosition either 'top','right','bottom' or 'left'. Where to position the legend.
#' @param out switch to chose the type of object returned by the function;
#'            default 'DF', i.e. the data.frame with the formatted data used for the plotting.
#'            Other options are:
#'                               'legend' for the legend of the plot,
#'                               'plot' for the figure itself,
#'                               'plot_noLegend' for the plot without the legend,
#'                               'plot_data' for the ggtable data to recreate the plot
#' @param print boolean; default: TRUE. Whether the plot should be printed
#'
#' @return depends on 'out' :
#'         'DF' returns the data.frame with all the required information for the plot,
#'         'legend' the legend of the plot,
#'         'plot' the figure itself,
#'         'plot_noLegend', the plot without the legend,
#'         'plot_data' the ggtable data to recreate the plot
#' @export
#'
nbGenesPerCond.unique<-function(DF, reorder=TRUE, xlab='Tissue',ylab,title,legend.title="Gene is",
                                publish=TRUE,threshold=0,thm,Y,remove.legend=FALSE,
                                legendPosition='bottom', out='DF',print=TRUE){

  uniqNoUniq=setNames(gg_color_hue(2),c('Tissue/cell specific','Unspecific to tissue or cell'))

  colNames<-colnames(DF) #Genes.ID are going to be added as a column
  names(colNames)<-colnames(DF) #so tempo would be a named list

  DF$Genes.ID<-rownames(DF)

  tempo<-parallel::mclapply(colNames,function (x){
    VEC<-DF[,x]
    names(VEC)<-rownames(DF)
    if(threshold==0){
      VEC<-VEC[VEC>0]
    }else{
      VEC<-VEC[VEC %>=% threshold]
    }
  })

  DF<-melt(DF,id.vars='Genes.ID')
  if(threshold==0){
    DF<-DF[DF$value>0,]
  }else{
    DF<-DF[DF$value %>=% threshold,]
  }
  DF$unique<-'Unspecific to tissue or cell'

  DF<-Reduce(rbind,lapply(colNames,function (x) {
    DIFF<-base::setdiff(colNames,x) #all the other conditions that will be tested against

    #next line gives the names that are unique to the current condition x
    VEC<-base::setdiff(names(tempo[[x]]),
                       base::unique(Reduce(c,lapply(DIFF,function(x) names(tempo[[x]])))) )
    if( dim(DF[DF$Genes.ID %in% VEC &DF$variable==x,])[1]>0){
      DF[DF$Genes.ID %in% VEC &DF$variable==x,]$unique<-'Tissue/cell specific'
    }
    DF[DF$variable==x,]
  }))

  DF$unique<-as.factor(DF$unique)
  if('Tissue/cell specific' %in% DF$unique)  DF$unique<-relevel(DF$unique,'Tissue/cell specific')

  if(reorder){
    v.order<-sapply(colNames,function(x) length(tempo[[x]]))
    #so the conditions are sorted by number of expressed genes.
    colNames<-names(v.order[base::order(v.order,decreasing=TRUE)])
    DF$variable<-factor(DF$variable,levels=colNames)#by
  }

  p <- ggplot(DF, aes_string(x="variable"))+geom_bar(aes(y=..count..,fill=unique))
  p <- p + scale_fill_manual(values=uniqNoUniq,drop=FALSE)
  if(publish) p <- p+theme_bw()
  if(!missing(xlab)) p <- p + xlab(xlab)
  if(!missing(ylab)) p <- p + ylab(ylab)
  if(!missing(Y)) p <- p + coord_cartesian(ylim=Y)
  if(!missing(thm)) p<-p+thm
  if(!missing(legend.title)) p <- p + guides(fill=guide_legend(title=legend.title))
  if(!missing(title)) p <- p + ggtitle(title)
  p <- p + theme(legend.position=legendPosition)

  p <- p +  theme(axis.text.x=element_text(angle=45,hjust = 1),plot.title=element_text(hjust=0.5))
  if(remove.legend) p <- p + theme(legend.position='none')

  if(print) print(p)

  switch(out,
         'DF'            = return(DF),
         'legend'        = return(g_legend(p)),
         'plot'          = return(print(p)),
         'plot_noLegend' = return(print(p+theme(legend.position='none'))),
         'plot_data'     = return(p))
}



#' Facetted plot (e.g. per study) of the distribution of genes across the columns.
#' Genes found uniquely in one condition are highlighted in another colour.
#'
#' @param DF numerical data.frame. Rownames have to be meaningful. Collates many columns from several studies.
#' @param xlab character string; used as x-axis label.
#' @param ylab character string; used as y-axis label.
#' @param title character string; title of the plot.
#' @param legend.title character string; used for the colour legend title
#' @param publish boolean. Default: TRUE.
#' @param threshold numeric; default: 0. Minimal level of expression to be considered as an expressed gene
#' @param thm other theme object (\code{\link[ggplot2]{theme}}) to pass on to the plot
#' @param Y numeric vector of two elements. Specify the limites of the y-axis.
#'          'c(minimal limit value, maximal limit value)'
#' @param remove.legend boolean. Default: FALSE.
#' @param out switch to chose the type of object returned by the function;
#'            default 'DF', i.e. the data.frame with the formatted data used for the plotting.
#'            Other options are:
#'                               'legend' for the legend of the plot,
#'                               'plot' for the figure itself,
#'                               'plot_noLegend' for the plot without the legend,
#'                               'plot_data' for the ggtable data to recreate the plot
#'                               'plot_def' for the output of \code{\link[ggplot2]{ggplot_build}} on the plot
#' @param print boolean; default: TRUE. Whether the plot should be printed
#' @param facetting character string; name of the column that should be used for the facetting.
#' @param flip boolean. Default= FALSE. Whether the x-axis and y-axis should be flipped or not.
#' @return depends on out. Objects/plot similar to \code{\link[barzinePhdR]{nbGenesPerCond.unique}}
#'                        with facetting for the different studies in the original dataframe.
#' @export
#'
complement2nbGenesPerCond.unique<-function(DF, xlab,ylab,title,legend.title,
                                           publish=TRUE,threshold=0,thm,Y,remove.legend=FALSE,
                                           out='DF',print=TRUE,facetting='Study',flip=FALSE){


  uniqNoUniq=setNames(gg_color_hue(2),c('Tissue/cell specific','Unspecific to tissue or cell'))

  p <- ggplot(DF, aes_string(x="variable")) + geom_bar(aes(y=..count..,fill=unique))
  p <- p + scale_fill_manual(values=uniqNoUniq,drop=FALSE)
  if(flip) {
    p <- p + theme_bw()+coord_flip()
    #p <- p + facet_grid(vars(Study),scales = 'free_y')
    p <- p + facet_grid(vars(as.formula(paste(facetting))),scales = 'free_y')
  }else{
    #p <- p + facet_grid(vars(Study))
    p <- p + facet_grid(vars(as.formula(paste(facetting))))
  }
  p <- p + theme(legend.position = 'bottom')

  if(publish) p <- p+theme_bw()
  if(!missing(xlab)) p <- p + xlab(xlab)
  if(!missing(ylab)) p <- p + ylab(ylab)
  if(!missing(Y)) p <- p + coord_cartesian(ylim=Y)
  if(!missing(thm)) p<-p+thm
  if(!missing(legend.title)) p <- p + guides(fill=guide_legend(title=legend.title))
  if(!missing(title)) p<- p + ggtitle(title)

  if(remove.legend) p <- p + theme(legend.position='none')
  p <- p + theme(axis.text.x=element_text(angle=35,hjust = 1),plot.title=element_text(hjust=0.5))

  if(print) print(p)

  switch(out,
         'DF'            = return(DF),
         'legend'        = return(g_legend(p)),
         'plot'          = return(print(p)),
         'plot_noLegend' = return(print(p+theme(legend.position='none'))),
         'plot_data'     = return(p),
         'plot_def'      = return(ggplot_build(p))
         )
}

#' Wrapper around gpplot2's density plot with preset arguments
#'
#' @param DF numeric data.frame containing expresion data
#' @param log2 boolean; default: TRUE.
#'             Whether a log2 transformation should be applied to the values
#' @param pseudocount numeric. default: 0.
#'                   Whether a pseudocount should be added before the log2 transformation is applied.
#' @param valueLab character string representing the unit in which the gene expression is measured;
#'                 added
#' @param baseFont string character. Name of the font to use for the plot.
#'                 Default: "Linux Libertine"
#' @param condCol character string vector. Colour palette for the plot.
#' The vector's names are for the conditions associated to each colour (content of the vector)
#' @param fontSize font size. Default: 11.
#' @param publish boolean. Default: TRUE. Whether to apply theme_bw() to the plot.
#' @param xmin numeric. Lower limit for the x-axis.
#' @param xmax numeric. Greater limitefor the x-axis.
#' @param ymin numeric. Lower limit for the y-axis.
#' @param ymax numeric. Greater limit for the y-axis.
#' @param removeLegend boolean, default: FALSE. Whethere to remove the legend
#' @param thm function theme( ) with the list of additional elements to apply
#' @param ylab character string. Label for the y-axis
#' @param alphaScale new alpha level in [0,1]. Default: 0.2
#'
#' @return a density plot
#' @export
#'
density_plot<-function(DF,log2=TRUE,pseudocount=0,valueLab='PSM',
                       baseFont='Linux Libertine',condCol,fontSize=11,
                       publish=TRUE,
                       xmin,xmax,ymin,ymax,removeLegend=FALSE,thm,
                       alphaScale=0.2,
                       ylab){
  DF_lg<-melt(DF)
  colnames(DF_lg)<-c('Tissue','value')

  if(log2) {
    if(pseudocount!=0)
      valueLab=paste0('Log2(',valueLab,' + ',pseudocount,')')
    else
      valueLab=paste0('Log2(',valueLab,')')
  }

  if(log2){
    p <- ggplot(DF_lg,aes(log2(value+pseudocount),fill=Tissue,colour=Tissue)) + geom_density(alpha=alphaScale)
  }else{
    p <- ggplot(DF_lg,aes(value,fill=Tissue,colour=Tissue)) + geom_density(alpha=alphaScale)
  }
  p <- p + xlab(valueLab)
  if(!missing(ylab)) p <- p + ylab(ylab)
  if(!missing(condCol)) p <- p + scale_fill_manual(values=condCol) + scale_colour_manual(values=condCol)
  p <- p + theme_bw(base_size = fontSize ,base_family = baseFont)
  if(!missing(xmin)&&!missing(xmax)) p <- p + xlim(xmin,xmax)
  if(!missing(ymin)&&!missing(ymax)) p <- p + ylim(ymin,ymax)
  if(removeLegend) p <- p + theme(legend.position='none')

  if(publish) p<- p + theme_bw()
  if(!missing(thm)) p<- p + thm
  p <- p+theme(plot.background = element_rect(fill = "transparent"))
  print(p)
}

#' Plot the distribution for the columns of a data.frame according to different geometries (ggplot2)
#'
#' @param DF numerical data.frame
#' @param transformation boolean; default: 'log2'. Whether if the data should be log2 transformed before being plotted
#' @param strip0 boolean; default: TRUE; Whether the 0 should be stripped from the data before plotting
#' @param geom character string; default 'all'. Allows to chose between the different available geometries.
#'             'all' output all available geometries. 'histogram' is for histogram, 'density' gives density plots.
#'             Anything else gives a boxplot.
#' @param all boolean; default: TRUE. Whether all geometries should be plotted.
#' @param debug boolean. Default: FALSE
#' @param verbose boolean. Default: FALSE
#' @param colourVal character vector. Colour palette.
#' @param limites numeric vector. Allows to define the minimal and maximal boundaries to apply to the different plots.
#' @param publish boolean. Whether to apply  \code{\link[ggplot2]{theme_bw}} to the plot(s).
#' @param valueLab Character vector; allows to provide custom units if needed
#'
#' @return a distribution plot of the input data.frame
#' @export
#'
plot_distribDF<-function(DF,transformation='log2',strip0=TRUE,geom='all',
                         all=FALSE,debug=FALSE,verbose=FALSE,colourVal,limites,
                         publish=TRUE,valueLab){

  DF_lg<-melt(DF)
  colnames(DF_lg)<-c('Tissue','value')

  if(debug){
    print(paste('nrow:',nrow(DF)))
    print(paste('str:',utils::str(DF_lg)))
  }

  if(all|!strip0){
    if(geom=='all'|all){
      if(missing(valueLab)) valueLab<-'expression'
      if(transformation=='log2'){
        if (verbose) print('log2 transformation applied')
        DF_lg$value<-log2(DF_lg$value)
        valueLab<-paste0('log2(',valueLab,')')
      }
      #histogram
      if (verbose)  print('Histogram with the 0')
      p <- ggplot(DF_lg, aes(x=value, y=..density..))
      p <- p + geom_histogram()+facet_wrap(~Tissue)+xlab(valueLab)
      if(publish) p <- p + theme_bw()
      p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
      if(!missing('limites')) p <- p + coord_cartesian(xlim=limites)
      print(p)

      #density
      if (verbose) print('Density plot with the 0')
      p <- ggplot(DF_lg,aes(value,fill=Tissue,colour=Tissue)) + geom_density(alpha=0.2)
      p <- p + xlab(valueLab)
      if(!missing(colourVal))
        p <- p + scale_fill_manual(values=colourVal) + scale_colour_manual(values=colourVal)

      if(publish) p <- p + theme_bw()
      if(!missing('limites')) p <-p + coord_cartesian(xlim=limites)
      print(p)

      #boxplot
      if (verbose) print('Boxplot with the 0')
      p <- ggplot(DF_lg, aes(Tissue,value),fill=Tissue)+ylab(valueLab)
      p <- p+geom_boxplot()
      if(!missing(colourVal))
        p <- p + scale_colour_manual(values=colourVal)+scale_fill_manual(values=colourVal)
      if(publish) p <- p + theme_bw()
      if(!missing('limites')) p <- p + coord_cartesian(ylim=limites)
      p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
      print(p)

    }else{
      if(geom=='histogram'){
        if(transformation=='log2'){
          if (verbose) print('log2 transformation applied')
          p <- ggplot(DF_lg, aes(x=log2(value), y=..density..))
        }else{
          if (verbose) print('No transformation applied')
          p <- ggplot(DF_lg, aes(x=log2(value+1), y=..density..))
        }
        if (verbose) print('histogram with the 0')
        p <- p + geom_histogram()+facet_wrap(~Tissue)
        if(publish) p<- p +theme_bw()
        p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
        if(!missing('limites')) p<-p + coord_cartesian(xlim=limites)
        print(p)
      }else{
        if(geom=='density'){
          if(transformation=='log2'){
            if (verbose) print('log2 transformation applied')
            p <- ggplot(DF_lg, aes(log2(value), fill = Tissue,colour=Tissue))
          }else{
            if (verbose) print('No transformation applied')
            p <- ggplot(DF_lg, aes(log2(value+1), fill = Tissue, colour=Tissue))
          }
          if (verbose) print('density with the 0')
          p <- p + geom_density(alpha = 0.2)
          if(publish) p<- p+ theme_bw()
          if(!missing('limites')) p<-p + coord_cartesian(xlim=limites)
          print(p)
        }else{
          #geom=='boxplot'
          if(transformation=='log2'){
            if (verbose) print('log2 transformation applied')
            p <- ggplot(DF_lg, aes(Tissue,log2(value)))
          }else{
            if (verbose) print('No transformation applied')
            p <- ggplot(DF_lg, aes(Tissue,value))
          }
          if (verbose) print('boxplot with the 0')
          p <- p+theme_bw()
          p <- p+geom_boxplot()+ theme(axis.text.x=element_text(angle=45,hjust=1))
          if(!missing('limites')) p<-p + coord_cartesian(ylim=limites)
          print(p)
        }
      }
    }
  }

  if(strip0 | all){
    if (verbose) print('all 0 stripped')
    DF_lg_n0<-DF_lg[DF_lg$value>0,]
    if(debug){
      print(paste('str:',utils::str(DF_lg)))
    }

    if(geom=='all'|all){
      if(missing(valueLab)) valueLab<-'expression'
      if(transformation=='log2'){
        if(verbose) print('log2 transformation applied')
        DF_lg_n0$value<-log2(DF_lg_n0$value)
        valueLab<-paste0('log2(',valueLab,')')
      }
      #histogram
      if (verbose) print('histogram without 0')
      p <- ggplot(DF_lg_n0, aes(x=value, y=..density..))
      p <- p + geom_histogram()+facet_wrap(~Tissue)+xlab(valueLab)
      if(publish) p <-p +theme_bw()
      p <- p + theme(axis.text.x=element_text(angle=45,hjust = 1))
      if(!missing('limites')) p<-p + coord_cartesian(xlim=limites)
      print(p)

      #density
      if (verbose) print('density without 0')
      p <- ggplot(DF_lg_n0,aes(value,fill=Tissue,colour=Tissue)) + geom_density(alpha=0.2)
      p <- p + xlab(valueLab)
      if(!missing(colourVal))
        p <- p + scale_fill_manual(values=colourVal) + scale_colour_manual(values=colourVal)

      if(publish) p<- p+theme_bw()
      if(!missing('limites')) p<-p + coord_cartesian(xlim=limites)
      print(p)

      #boxplot
      if (verbose) print('boxplot without 0')
      p <- ggplot(DF_lg_n0, aes(Tissue,value))+ylab(valueLab)
      p <- p+geom_boxplot()
      if(!missing(colourVal))
        p <- p + scale_colour_manual(values=colourVal)+scale_fill_manual(values=colourVal)
      if(publish) p<- p +theme_bw()
      p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
      if(!missing('limites')) p <-p + coord_cartesian(ylim=limites)
      print(p)
    }else{
      if(geom=='histogram'){
        if(transformation=='log2'){
          if (verbose) print('log2 transformation applied')
          p<- ggplot(DF_lg_n0, aes(x=log2(value), y=..density..))
        }else{
          if (verbose) print ('no transformation applied')
          p<- ggplot(DF_lg_n0, aes(x=log2(value+1), y=..density..))
        }
        if(verbose) print('histogram without 0')
        p <- p + geom_histogram()+facet_wrap(~Tissue)
        if(publish) p <- p +theme_bw()
        p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
        if(!missing('limites')) p <-p + coord_cartesian(xlim=limites)
        print(p)
      }else{
        if(geom=='density'){
          if(verbose) print('density without 0')
          if(transformation=='log2'){
            if(verbose) print('log2 transformation applied')
            p <- ggplot(DF_lg_n0, aes(log2(value), fill = Tissue, colour=Tissue))
          }else{
            if (verbose) print('no transformation applied')
            p <- ggplot(DF_lg_n0, aes(log2(value+1), fill = Tissue, colour=Tissue))
          }
          p <- p + geom_density(alpha = 0.2)
          if(publish) p <- p+theme_bw()
          if(!missing('limites')) p<-p +coord_cartesian(xlim=limites)
          print(p)
        }else{
          #geom=='boxplot'
          if(verbose) print ('boxplot without 0')
          if(transformation=='log2'){
            if (verbose) print('log2 transformation applied')
            p <- ggplot(DF_lg_n0, aes(Tissue,log2(value)))
          }else{
            if(verbose) print('no transormation applied')
            p <- ggplot(DF_lg_n0, aes(Tissue,value))
          }
          p <- p+geom_boxplot()
          if (publish) p<- p +theme_bw()
          p <- p + theme(axis.text.x=element_text(angle=45,hjust=1))
          if(!missing('limites')) p<-p + coord_cartesian(ylim=limites)
          print(p)
        }
      }
    }
  }
}


# Gene oriented -------------


#' Extract the genes that are specific to a tissue
#'
#' @param DF numeric data.frame
#' @param threshold numeric. Minimal expression to be considered.
#'              For 0, the comparison is strict.
#' @param strip boolean. Default: TRUE.
#' @param format character string. 'lg' allows to output a long-format data.frame
#' @param rowNames boolean. Default: TRUE.
#' @param verbose boolean. Default: FALSE.
#'
#' @return a long or wide data.frame
#' @export
#'
extractSpe<-function(DF,threshold,strip=TRUE,format='lg',rowNames=TRUE,
                     verbose=FALSE){

  if(strip|rowNames) {
    if(format!='lg') {
      print("Due to other selected arguments, format switched to long ('lg')")
      format='lg'
    }}
  DF<-computeBreadth(DF,threshold=threshold,omit.zero=TRUE,typeR='df')
  DF.spe<-na.omit(DF[DF$nb.tissues==1,!colnames(DF) %in% 'nb.tissues'])
  if(format=='lg') {
    if(rowNames) DF.spe$geneID<-rownames(DF.spe)
    DF.spe<-suppressMessages(melt(DF.spe))
    if(strip) {
      if(threshold==0){
        DF.spe<-na.omit(DF.spe[DF.spe$value > 0,])
      }else{
        DF.spe<-na.omit(DF.spe[DF.spe$value %>=% threshold,])
      }
    }
  }
  if(verbose)
    print(paste0(nrow(DF.spe),' features have been found in only once at ',threshold,' (threshold)'))
  return(DF.spe)
}


#' Compute breadth of expression of the genes across a data.frame
#'
#' @param DF numeric data.frame
#' @param threshold numeric. Minimal expression to be considered.
#'              For 0, the comparison is strict.
#' @param omit.zero boolean.
#'              Whether the genes which do not reach the minimal expression threshold in any tissue
#'              should be stripped from the output.
#' @param typeR character string. Chose the type of result to return.
#'              'vec' for a vector of the observed expression breadths,
#'              'nameVec' for a named vector of the observed expression breadths
#'              associated to each rownames,
#'              'df' or 'DF' for the original data.frame with an additional column,
#'              'dfBool' for a boolean data.frame where TRUE when the gene expression reaches
#'              the minimal expression.
#'              'unique_tissueList' data.frame with genes that are found in one tissue only
#'               above the minimal threshold.
#'              'plot' for a histogram of the gene expression breadth in that data.frame
#'
#' @return the expression breadth of the genes in different (see typeR)
#' @export
#'
computeBreadth<-function(DF,threshold,omit.zero,typeR='vec'){
  if(typeR!='dfBool'){
    if(threshold==0){
      DF$nb.tissues<-rowSums(DF > threshold)
    }else{
      DF$nb.tissues<-rowSums(DF %>=% threshold)
    }
    if(typeR=='unique_tissueList'){
      #only genes expressed in one tissue
      #and new column that gives the name of the tissue
      DF<-na.omit(DF[DF$nb.tissues==1,])
      DF$nb.tissues<-NULL
      DF$condition<-as.character(lapply(row.names(DF),function(x){
        return(names(which.max(DF[x,])))
      }
      ))
    }
  }else{ #for paired tissues analysis in particular
    if(threshold==0){
      DF<-as.data.frame(DF>threshold)
    }else{
      DF<-as.data.frame(DF %>=% threshold)
    }
    DF$nb.tissues<-rowSums(DF)
  }
  if(omit.zero & typeR=='vec') warning('Even if the omit.zero option is selected, all the rows are returned')
  if(omit.zero & typeR=='df') DF<-na.omit(DF[DF$nb.tissues!=0,])
  if(omit.zero & typeR=='plot') DF<-na.omit(DF[DF$nb.tissues!=0,])
  switch(typeR,
         'vec'= return(as.numeric(DF$nb.tissues)),
         'nameVec'= {nameVec=setNames(DF$nb.tissues,rownames(DF))
         return(nameVec)
         },
         'df' = return(DF),
         'DF' = return(DF),
         'dfBool'= return(DF),
         'unique_tissueList'= return(DF),
         'plot'={
           p<-ggplot(DF,aes(as.factor(nb.tissues)))+geom_bar()
           p<-p + xlab('Number of tissues')+ ylab('Genes count')
           return(p)
         }
  )
}


#' More customisable plot based on computeBreadth
#'
#' @param DF numeric data.frame
#' @param omit.zero numeric. Minimal expression to be considered.
#'        For 0, the comparison is strict.
#' @param threshold numeric. Minimal expression to be considered.
#'        For 0, the comparison is strict.
#' @param annot boolean. Default: TRUE.
#' @param publish boolean. Default: TRUE.
#'        Whether to apply \code{\link[ggplot2]{theme_bw}} to the plot.
#' @param ... allow to handle more theme elements in \code{\link[ggplot2]{theme_bw}}
#'
#' @return a plot
#' @export
#'
plotUniq.distrib<-function(DF,omit.zero=TRUE,threshold=0,
                           annot=TRUE,
                           publish=TRUE,...){

  DF<-computeBreadth(DF,omit.zero=omit.zero,threshold=threshold,typeR='df')

  p<-ggplot(DF,aes(x=factor(nb.tissues)))+geom_bar()+labs(x='Number of tissues')
  if(annot)    p<-p+ stat_count(geom="text", aes(label=..count..),
                                vjust=-0.5, position="identity",size=4)
  if (publish) p<- p +theme_bw(...)
  return(p)
}

