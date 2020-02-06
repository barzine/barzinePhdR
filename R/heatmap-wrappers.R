#' Draw a customised heatmap with marked clusters in the margins for a correlation matrix
#'
#' @param data correlation data.frame to plot
#' @param label.A character string. Name of the first data.frame used for the correlation matrix
#'               for the additional text.
#' @param label.B character string. Name of the second data.frame used for the correlation matrix
#'              for the additional text.
#' @param out.A character string for the plot (x-axis).
#'              Name of the first data.frame used for the correlation matrix
#' @param out.B character string for the plot (y-axis).
#'              Name of the second data.frame used for the correlation matrix
#' @param col color palette to use; default barzinePhdR::colC(),
#'            i.e. grDevices::colorRampPalette('aliceblue','darkcyan')
#' @param dend character string that correspond to (gplots) heatmap.2 "dendrogram"
#'             indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms.
#'             Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both',
#'             then a warning is issued and Rowv (or Colv) arguments are honoured.
#' @param margins numeric vector of length 2 containing the margins
#'                (see par(mar= *)) for column and row names, respectively.
#' @param ... other parameters that can be handled by (gplots) heatmap.2
#' @param report boolean. Default: TRUE. Output additional text (example for an (r)html report)
#' @param index character string. default:"####" Allows to customise the level of indexing through markdown tags.
#'
#' @return a heatmap plot
#' @export
#'
heatmapWithClustering<-function(data,label.A,label.B,out.A,out.B,dend='none',
                                margins,col=barzinePhdR::colCbc(),
                                ...,report=TRUE,index="####"){

  if(missing(margins)) margins=c(10,10)
  if (report)
    cat(paste0("\n\n",index," Heatmap (with clustering) for correlation matrix of ",
               label.A," and ",label.B,":\n\n"))

  heatmap.2(t(as.matrix(data)),hclustfun=function(x) hclust(x,method = 'ward.D'),
            distfun=function(c) as.dist(1 - c), trace="none", dendrogram=dend,
            col=col,cexRow =1.4,cexCol = 1.4,margins=margins,
            xlab=out.A,ylab=out.B,...)
}



#' Draw a customised heatmap where the different variables are sorted alphabetically
#'
#' @param data correlation (numeric) data.frame to plot
#' @param label.A character string. Name of the first data.frame used for the correlation matrix
#'               for the additional text.
#' @param label.B character string. Name of the second data.frame used for the correlation matrix
#'              for the additional text.
#' @param out.A character string for the plot (x-axis).
#'              Name of the first data.frame used for the correlation matrix
#' @param out.B character string for the plot (y-axis).
#'              Name of the second data.frame used for the correlation matrix
#' @param digits default: 1. How many digits should be kept for the annotation.
#' @param col color palette to use; default barzinePhdR::colC(),
#'            i.e. grDevices::colorRampPalette('aliceblue','darkcyan')
#' @param key boolean. default: TRUE. Whether the key of the heatmap should also be drawn.
#' @param margins numeric vector of length 2 containing the margins
#'                (see par(mar= *)) for column and row names, respectively.
#' @param ... other parameters that can be handled by (gplots) heatmap.2
#' @param report boolean. Default: TRUE. Output additional text (example for an (r)html report)
#' @param index character string. default:"####" Allows to customise the level of indexing through markdown tags.
#' @param annot boolean; default: TRUE.
#'              Whether the input data.frame should be rordered alphabetically before being plotted.
#' @param sorting boolean; default: TRUE.
#'                Whether the input data.frame should be rordered alphabetically before being plotted.
#' @param notecol character string specifying the color for cellnote text. Defaults to "grey90".
#' @param srtCol angle of column labels, in degrees from horizontal. Default: 45
#' @param cexRow positive numbers, used as cex.axis in for the row axis labeling.
#'               Default: 1.4
#' @param cexCol positive numbers, used as cex.axis in for the column axis labeling.
#'               Default: 1.4
#' @param cex.axis positive numbers, used as cex.axis in for the row and column axis labeling.
#'               Default: 1.4
#' @param mathMethod character string that allows to pick how the correlation (in cellnote) should be rounded.
#'                   Default: "signif"
#'
#' @return a heatmap plot
#' @export
#'
heatmapAlphabet<-function(data,label.A,label.B,out.A,out.B,col=barzinePhdR::colCbc(),
                          digits=1,key,margins,...,report=TRUE,index="####",annot='diag',
                          sorting=TRUE,notecol="grey90",srtCol=45,
                          cexRow=1.4,cexCol=1.4, cex.axis=1.4,
                          mathMethod='signif'){

  if(missing(key)) key=TRUE
  if(missing(margins)) margins=c(10,10)
  if(missing(out.A)) out.A=label.A
  if(missing(out.B)) out.B=label.B
  if (report)
    cat(paste0("\n\n",index," Heatmap (no clustering, ordered by name) for correlation matrix of ",
               label.A," and ",label.B,":\n\n"))
  if(sorting){
    data<-data[order(rownames(data)),]
    data<-data[,order(colnames(data))]
  }

  if(mathMethod %in% c('signif','round')) {
    signif2<-function(x){
      signif(x,digits = digits)
    }
    round2<-function(x){
      round(x,digits= digits)
    }
  }

  cellAnnot<-apply(eval(call(name=paste0(mathMethod,'2'),as.matrix(data))),2,as.character)
  if(annot=='diag'){
    tmp<-diag(cellAnnot)
    cellAnnot<-matrix(data="",nrow=nrow(data),ncol(data))
    diag(cellAnnot)<-tmp

    heatmap.2(t(as.matrix(data)), trace="none", dendrogram="none",
              Rowv=FALSE,Colv=FALSE, col=col,
              margins=margins,xlab=out.A,ylab=out.B,key=key,
              cexRow=cexRow,cexCol=cexCol, cex.axis=cex.axis,
              cellnote=cellAnnot,notecol=notecol,srtCol=srtCol,...)
  }else{
    if(annot=='complete'){
      heatmap.2(t(as.matrix(data)), trace="none", dendrogram="none",
                Rowv=FALSE,Colv=FALSE, col=col,
                cexRow=cexRow,cexCol=cexCol, cex.axis=cex.axis,
                margins=margins,xlab=out.A,ylab=out.B,key=key,
                cellnote = cellAnnot,
                notecol=notecol,srtCol=srtCol,...)
    }else{
      heatmap.2(t(as.matrix(data)), trace="none", dendrogram="none",
                Rowv=FALSE,Colv=FALSE, col=col,
                cexRow=cexRow,cexCol=cexCol, cex.axis=cex.axis,
                margins=margins,xlab=out.A,ylab=out.B,srtCol=srtCol,...)
    }
  }

}


#' Draw a customised heatmap where the different variables are sorted alphabetically
#'
#' @param data correlation (numeric) data.frame to plot
#' @param label.A character string. Name of the first data.frame used for the correlation matrix
#'               for the additional text.
#' @param label.B character string. Name of the second data.frame used for the correlation matrix
#'              for the additional text.
#' @param out.A character string for the plot (x-axis).
#'              Name of the first data.frame used for the correlation matrix
#' @param out.B character string for the plot (y-axis).
#'              Name of the second data.frame used for the correlation matrix
#' @param digits integer. default: 1. How many digits should be kept for the annotation.
#' @param key boolean. default: TRUE. Whether the key of the heatmap should also be drawn.
#' @param col color palette to use; default barzinePhdR::colC(),
#'            i.e. grDevices::colorRampPalette('aliceblue','darkcyan')
#' @param margins numeric vector of length 2 containing the margins
#'                (see par(mar= *)) for column and row names, respectively.
#' @param ... other parameters that can be handled by (gplots) heatmap.2
#' @param report boolean. Default: TRUE. Output additional text (example for an (r)html report)
#' @param index character string. default:"####" Allows to customise the level of indexing through markdown tags.
#' @param annot character string. Default: 'diag'. Annotate only the diagonal of the heatmap
#'              with the value of correlation.
#' @param sorting boolean; default: TRUE.
#'                Whether the input data.frame should be rordered alphabetically before being plotted.
#'
#' @return a heatmap plot
#' @export
#'
heatmapAlpabet<-function(data,label.A,label.B,out.A,out.B,digits,key,
                         col=barzinePhdR::colCbc(),
                         sorting=TRUE,
                         margins,...,report=TRUE,index="####",annot='diag'){

  if(missing(digits)) digits=1
  if(missing(key)) key=TRUE
  if(missing(margins)) margins=c(10,10)
  if(missing(out.A)) out.A=label.A
  if(missing(out.B)) out.B=label.B

  if (report)
    cat(paste0("\n\n",index," Heatmap (no clustering, ordered by name) for correlation matrix of ",
               label.A," and ",label.B,":\n\n"))
  if(sorting){
    data<-data[order(rownames(data)),]
    data<-data[,order(colnames(data))]
  }

  cellAnnot<-apply(signif(as.matrix(data),digits=digits),2,as.character)


  if(annot=='diag'){
    tmp<-diag(cellAnnot)
    cellAnnot<-matrix(data="",nrow=nrow(data),ncol(data))
    diag(cellAnnot)<-tmp

    heatmap.2(t(as.matrix(data)), trace="none", dendrogram="none",
              Rowv=FALSE,Colv=FALSE, col=col,
              cexRow=1,cexCol=1, cex.axis=1, margins=margins,
              xlab=out.B,ylab=out.A,key=key,
              cellnote=cellAnnot,notecol="grey90",...)
  }else{
    heatmap.2(t(as.matrix(data)), trace="none", dendrogram="none",
              Rowv=FALSE,Colv=FALSE, col=col,
              cexRow =1.4,cexCol = 1.4,margins=margins,xlab=out.B,ylab=out.A,...)
  }

}

