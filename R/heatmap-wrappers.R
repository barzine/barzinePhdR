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



#' Allows to display the correct colour annotation
#'
#' @param Names vector of character strings. name of the columns
#' @param cond vector of character strings. Name of the tissues/conditions to map the columns
#' @param palette named vector palette. (names are the same as condition)
#'
#' @return a properly formatted vector of colours
#' @export
#'
matchColCond<-function(Names,cond,palette){
  stopifnot(all(cond %in% names(palette)))
  res<-lapply(Names,function(x){
    temp<-lapply(cond,function(y){
      ifelse(grep(y,x,ignore.case = TRUE),return(palette[y]),return(NA))
    })
    temp<-temp[!is.na(temp)]
    return(temp)
  })
  return(unlist(res))
}

#' Wrapper around heatmap.2; allows to annotate with the colours on the side.
#'
#' @param DF data.frame or matrix to use for the heatmap
#' @param method string to pick the method with which to compute the correlation.
#'               One of "pearson" (default), "kendall", or "spearman"
#' @param use character string giving a method for computing covariances in the presence of missing values.
#'            This must be (an abbreviation of) one of the strings
#'            "everything" (default), "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param hclustMethod character string. he agglomeration method to be used.
#'                     This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param trace character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'none'.
#' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured.
#' @param col colors used for the image. Defaults to col=colorRampPalette(c('ghostwhite','darkcyan'))
#' @param cexRow,cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
#' @param annotate boolean. Whether to add the correlations as annotation on the heatmap
#' @param signifdigit integer. Number of digit to show for the correlation showed as annotation. Default:1
#' @param rowsidecolors (optional) character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
#'                       Otherwise can be created internally
#' @param colsidecolors (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#'                      Otherwise can be created internally
#' @param density.info character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key.
#' @param notecol character string. Colour for the annotation
#' @param ColsideMode character string that allows to create the annotation colours based on the column names
#'                    one of the following strings: 'match', 'extract' or 'extractSpace'
#' @param RowsideMode character string that allows to create the annotation colours based on the row names (of the correlation matrix)
#'                    one of the following strings: 'extract', 'extractSpace', 'extract2' or 'extractSpace2'
#' @param srtCol angle of column labels, in degrees from horizontal. Default: 45
#' @param margins numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively. Default: c(8,12)
#' @param common.cond argument for matchColCond
#' @param datasetCol colour palette for the datasets in the form of a named vector
#' @param TissueCol colour palette for the tissues in the form of a named vector
#' @param baseFont character string. Allows to pick the font of the figure
#' @param ... other parameters that can be used by heatmap.2
#' @param parenthesis boolean. Whether the column names comprises parenthesis
#' @param key.title main title of the color key. If set to NA no title will be plotted.
#'
#' @return a heatmap
#' @export
#'
mHeatmap<-function(DF,method='pearson',use='everything',hclustMethod='ward.D',
                  trace='none',dendrogram='row',col,cexRow=1,cexCol=1,
                  annotate=FALSE,signifdigit=1, rowsidecolors,colsidecolors,density.info,
                  notecol='black',ColsideMode='match',RowsideMode='extract',srtCol=45,
                  margins=c(12,8),common.cond,datasetCol,TissueCol,baseFont,...,parenthesis=TRUE,key.title){

  if(!missing(baseFont)) op <- par(family = baseFont)

  if (parenthesis) names(datasetCol)<-sapply(names(datasetCol), function(x) paste0('(',x,')'))
  if(missing(key.title)) key.title=NULL

  if(missing(col)) col=grDevices::colorRampPalette(c('ghostwhite','darkcyan'))

  corDF<-cor(DF,method=method,use=use)

  if(missing(rowsidecolors)){
    if(RowsideMode=='extract')
      rowsidecolors=datasetCol[sapply(rownames(corDF),
                                      function(x){return(unlist(strsplit(x,'\\.')[[1]][2]))})]
    if(RowsideMode=='extractSpace')
      rowsidecolors=datasetCol[sapply(rownames(corDF),
                                      function(x){return(unlist(strsplit(x,' ')[[1]][2]))})]
    if(RowsideMode=='extract2')
      rowsidecolors=datasetCol[sapply(rownames(corDF),function(x){
        return(unlist(strsplit(x,'\\.')[[1]][length(unlist(strsplit(x,'\\.')))]))})]
    if(RowsideMode=='extractSpace2')
      rowsidecolors=datasetCol[sapply(rownames(corDF),function(x){
        return(unlist(strsplit(x,' ')[[1]][length(unlist(strsplit(x,' ')))]))})]
  }

  if(missing(colsidecolors)){
    if(ColsideMode=='match')
      colsidecolors=matchColCond(colnames(corDF),common.cond,TissueCol)
    if(ColsideMode=='extract')
      colsidecolors=TissueCol[sapply(rownames(corDF),
                                     function(x){
                                       return(unlist(strsplit(x,'\\.')[[1]][1]))})]
    if(ColsideMode=='extractSpace')
      colsidecolors=TissueCol[sapply(rownames(corDF),
                                     function(x){
                                       return(unlist(strsplit(x,' ')[[1]][1]))})]
  }

  if(annotate){
    heatmap.2(corDF,
              hclustfun = function(x) hclust(x,method=hclustMethod),
              distfun = function(c) as.dist(1 - c), trace=trace,
              dendrogram = dendrogram, col=col,
              cexRow = cexRow, cexCol = cexCol,
              cellnote=signif(corDF,signifdigit),notecol = 'black',
              RowSideColors = rowsidecolors,
              ColSideColors = colsidecolors,
              srtCol=srtCol,margins=margins,key.xlab = key.title, key.title=NA)

  }else{
    heatmap.2(corDF,
              hclustfun = function(x) hclust(x,method=hclustMethod),
              distfun = function(c) as.dist(1 - c), trace=trace,
              dendrogram = dendrogram, col=col,
              cexRow = cexRow, cexCol = cexCol,
              RowSideColors = rowsidecolors,
              ColSideColors = colsidecolors,
              srtCol=srtCol,margins=margins,key.xlab = key.title, key.title=NA,...)
  }
}


#' Wrapper around heatmap.2; allows to annotate with the colours on the side.
#' Colours are the transpose of mHeatmap
#'
#' @param DF data.frame or matrix to use for the heatmap
#' @param method string to pick the method with which to compute the correlation.
#'               One of "pearson" (default), "kendall", or "spearman"
#' @param use character string giving a method for computing covariances in the presence of missing values.
#'            This must be (an abbreviation of) one of the strings
#'            "everything" (default), "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param hclustMethod character string. he agglomeration method to be used.
#'                     This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single", "complete", "average" (= UPGMA), "mcquitty" (= WPGMA), "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param trace character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'none'.
#' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured.
#' @param col  colors used for the image. Defaults to col=colorRampPalette(c('ghostwhite','darkcyan'))
#' @param cexRow,cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
#' @param annotate boolean. Whether to add the correlations as annotation on the heatmap
#' @param signifdigit integer. Number of digit to show for the correlation showed as annotation. Default:1
#' @param rowsidecolors (optional) character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
#'                       Otherwise can be created internally
#' @param colsidecolors (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#'                      Otherwise can be created internally
#' @param notecol character string. Colour for the annotation
#' @param ColsideMode character string that allows to create the annotation colours based on the column names
#'                    one of the following strings: 'extract', 'extractSpace', 'extract2' or 'extractSpace2'
#' @param RowsideMode character string that allows to create the annotation colours based on the row names (of the correlation matrix)
#'                    one of the following strings: 'match', 'extract' or 'extractSpace'
#' @param srtCol angle of column labels, in degrees from horizontal. Default: 45
#' @param margins numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively. Default: c(8,12)
#' @param common.cond argument for matchColCond
#' @param datasetCol colour palette for the datasets in the form of a named vector
#' @param TissueCol colour palette for the tissues in the form of a named vector
#' @param baseFont character string. Allows to pick the font of the figure
#' @param ... other parameters that can be used by heatmap.2
#' @param parenthesis boolean. Whether the column names comprises parenthesis
#' @param key.title main title of the color key. If set to NA no title will be plotted.
#'
#' @return a heatmap
#' @export
#'
mHeatmapTC<-function(DF,method='pearson',use='pairwise.complete.obs',hclustMethod='ward.D',
                    trace='none',dendrogram='row',col,cexRow=1,cexCol=1,
                    annotate=FALSE,signifdigit=1, rowsidecolors,colsidecolors,
                    notecol='black',ColsideMode='match',RowsideMode='extract',srtCol=45,
                    margins=c(12,8),common.cond,datasetCol,TissueCol,baseFont,...,parenthesis=TRUE,key.title){

  if(!missing(baseFont)) op <- par(family = baseFont)
  if (parenthesis) names(datasetCol)<-sapply(names(datasetCol), function(x) paste0('(',x,')'))
  if(missing(key.title)) key.title=NULL

  if(missing(col)) col=grDevices::colorRampPalette(c('ghostwhite','darkcyan'))

  corDF<-cor(DF,method=method,use=use)

  if(missing(colsidecolors)){
    if(ColsideMode=='extract')
      colsidecolors=datasetCol[sapply(rownames(corDF),
                                      function(x){return(unlist(strsplit(x,'\\.')[[1]][2]))})]
    if(ColsideMode=='extractSpace')
      colsidecolors=datasetCol[sapply(rownames(corDF),
                                      function(x){return(unlist(strsplit(x,' ')[[1]][2]))})]
    if(ColsideMode=='extract2')
      colsidecolors=datasetCol[sapply(rownames(corDF),function(x){
        return(unlist(strsplit(x,'\\.')[[1]][length(unlist(strsplit(x,'\\.')))]))})]
    if(ColsideMode=='extractSpace2')
      colsidecolors=datasetCol[sapply(rownames(corDF),function(x){
        return(unlist(strsplit(x,' ')[[1]][length(unlist(strsplit(x,' ')))]))})]
  }

  if(missing(rowsidecolors)){
    if(RowsideMode=='match')
      rowsidecolors=matchColCond(colnames(corDF),common.cond,TissueCol)
    if(RowsideMode=='extract')
      rowsidecolors=TissueCol[sapply(rownames(corDF),
                                     function(x){
                                       return(unlist(strsplit(x,'\\.')[[1]][1]))})]
    if(RowsideMode=='extractSpace')
      rowsidecolors=TissueCol[sapply(rownames(corDF),
                                     function(x){
                                       return(unlist(strsplit(x,' ')[[1]][1]))})]
  }

  if(annotate){
    heatmap.2(corDF,
              hclustfun = function(x) hclust(x,method=hclustMethod),
              distfun = function(c) as.dist(1 - c), trace=trace,
              dendrogram = dendrogram, col=col,
              cexRow = cexRow, cexCol = cexCol,
              cellnote=signif(corDF,signifdigit),notecol = 'black',
              RowSideColors = rowsidecolors,
              ColSideColors = colsidecolors,
              srtCol=srtCol,margins=margins,key.xlab = key.title, key.title=NA, ...)

  }else{
    heatmap.2(corDF,
              hclustfun = function(x) hclust(x,method=hclustMethod),
              distfun = function(c) as.dist(1 - c), trace=trace,
              dendrogram = dendrogram, col=col,
              cexRow = cexRow, cexCol = cexCol,
              RowSideColors = rowsidecolors,
              ColSideColors = colsidecolors,
              srtCol=srtCol,margins=margins,key.xlab = key.title, key.title=NA, ...)
  }
}
