## Analyses based on correlation ------------------------------------------------

#' Compute the correlation between two matrices variables( alphabetical order)
#'
#' @param A numeric data.frame for the first study to compare
#' @param B numeric data.frame for the second study to compare
#' @param label.A character string. Name of the first study used in the report.
#' @param label.B character string. Name of the second study used in the report.
#' @param method character string. Method for base::cor, either "pearson", "spearman" or
#'              "kendall".
#' @param log2 logical. Default: TRUE. Whether the data should be transform with log2 or not
#' @param report boolean. Default:TRUE. Whether to output the result before returning the matrix.
#' @param Latex boolean. Default: FALSE. Whether to output in latex format.
#' @param Kable boolean. Default: TRUE. Whether to report with kable.
#' @param Grid boolean. Default: TRUE. Whether to create an image of the output.
#' @param headers boolean. Default: TRUE. Whether to use headers for the reporting.
#'
#' @return a matrix with the columns correlation of the two studies.
#' @export
#'
compute_tissue_corr_matrix<-function(A,B,label.A,label.B,method,log2=TRUE,
                                     report=TRUE,Latex=FALSE,Kable=TRUE,Grid=TRUE,
                                     headers=TRUE){

  if (report){
    if(headers)  cat(paste0("\n\n### Correlation between ",label.A," and ",label.B,":\n\n"))
  }
  A<-A[,sort(colnames(A))]
  B<-B[,sort(colnames(B))]
  if(log2){
  CorrMat<-data.frame(lapply(colnames(B),function(x){
    sapply(colnames(A),function(y){
      cor(log2(A[,y]+1),log2(B[,x]+1),method=method)})}))
  }else{
    CorrMat<-data.frame(lapply(colnames(B),function(x){
      sapply(colnames(A),function(y){
        cor(A[,y],B[,x],method=method)})}))
  }

  colnames(CorrMat)<-colnames(B)
  if(report){
    if(Latex) xtable::xtable(CorrMat)
    if(Kable) print(knitr::kable(CorrMat))
    if(Grid){
      gridExtra::grid.table(signif(CorrMat,digits=2))
      grid::grid.newpage()
    }
  }
  return(CorrMat)
}

#' Create the data.frame and can also test and plot
#' the comparison between the tissues of the first study with the second one
#'
#' @param CorrMat correlation matrix between two studies tissue expression (column correlation)
#' @param test boolean. Default: FALSE. Whether a (Welch's) t-test should be performed
#'             between the group where the same tissue from each study are compared together
#'             and the group where one tissue from the first study is compared to
#'             a random tissue of the second tissue.
#' @param plot boolean. Default: FALSE. Whether the \strong{Test} should be plotted as well.
#' @param publish boolean. Default:TRUE. Whether the theme theme_bw() should be applied to the plot
#' @param x numeric. Default:0.75. x-axis of where to center the p-value.
#' @param y numeric. Default:0.1. y-axis of where to center the p-value.
#' @param ym numeric. Default: 0.02. y-axis of where to write the medium value of the correlation
#' @param textaxis numeric. Default:15. Size of the text on the x-axis
#' @param titleaxis  numeric. Default:15. Size of the x-axis' title.
#' @param out character string. Allows to pick the output.
#'            "DFtest" (default):
#' @param sep character or character string. Default: '.'
#' @param remove.identity boolean. Default:TRUE.
#'        Whether the cases where the same tissue from the first study and the second one should be removed.
#' @param reorder boolean. Default: TRUE. Whether the data.frame should be order in decreasing order of the observed correlation
#'
#' @return depends on 'out'. Can be a data.frame, a plot, only the identical tissue comparisons,
#'         only the different tissue comparisons, ...
#' @export
#'
test_CorrMat_identical_different<-function(CorrMat,test=FALSE,plot=FALSE,publish=TRUE,
                                           x=0.75,y=0.1,ym=0.02,textaxis=15,titleaxis=15,
                                           out='DFtest',sep='.',remove.identity=TRUE,reorder=TRUE){

  DFtest<-as.data.frame(t(data.frame(utils::combn(colnames(CorrMat),2,simplify=FALSE))),stringsAsFactors = FALSE)
  rownames(DFtest)<-1:nrow(DFtest)
  colnames(DFtest)<-c('Pair1','Pair2')

  sep=paste0('[',sep,']')
  DFtest$Type<-sapply(1:nrow(DFtest),
                      function(x) ifelse(strsplit(DFtest$Pair1[x],sep)[[1]][1]==strsplit(DFtest$Pair2[x],sep)[[1]][1],
                                         'Same tissues','Different tissues'))
  DFtest$Correlation<-sapply(1:nrow(DFtest),function(x) CorrMat[DFtest[x,1],DFtest[x,2]])
  if(remove.identity) DFtest<-DFtest[DFtest$Pair1!=DFtest$Pair2,]
  if(reorder) DFtest<-DFtest[base::order(DFtest$Correlation,decreasing=TRUE),]

  if(test){
    Identical<-DFtest[DFtest$Type=="Same tissues",]
    Different<-DFtest[DFtest$Type=="Different tissues",]
    test<-t.test(Identical$Correlation,Different$Correlation)
    if(reorder) Identical<-Identical[base::order(Identical$Correlation,decreasing=TRUE),]
    estimateDF<-data.frame(Type=c("Same tissues","Different tissues"),
                           Correlation=test$estimate)
    if(plot){
      p <- ggplot2::ggplot(DFtest,aes(Type,Correlation))+ggplot2::geom_boxplot()+ggplot2::coord_cartesian(ylim=c(0,1))
      p <- p + ggplot2::annotate('text',label=paste('p-value=',signif(test$p.value,digits=8)),
                                 x=x,y=y)+ggplot2::geom_text(data=estimateDF,
                                                             aes(label=paste0('m = ',
                                                                              signif(Correlation,digits=4)),
                                                                 y=ym),
                                                             size=10)
      if(publish) p <- p + theme_bw()
      p <- p + ggplot2::theme(axis.text=element_text(size=textaxis),
                              axis.title=element_text(size=titleaxis),
                              legend.position='none')

      print(p)
      switch(out,
             'identical'= return(Identical),
             'different'= return(Different),
             'DFtest'= return(DFtest),
             'test' = return(test),
             'plot' = return(p),
             'test_identical'= return(list(test,identical)),
             'DFtest_identical'= return(list(DFtest,identical))
      )
    }
  }else{
    return(DFtest)
  }
}


#' Test with a Welch t-test between the group of similar tissues versus the group of different tissues.
#'
#' @param CorrMat correlation matrix between two studies tissue expression (column correlation)
#' @param label.A character string. Name of the first study on which the correlations are based
#' @param label.B character string. Name of the second study on which the correlations are based
#' @param report boolean. Default: TRUE. Whether to output a text report.
#' @param index character string. For the notebook, allows to create new sections.
#' @param output character string.
#'               'identical' to return only the pairs that comprise the same tissues,
#'               'different' to return only the pairs that comprise the different tissues,
#'               'DFtest' to return all the pairs,
#'               'test' returns only the test,
#'               'plot' returns only the plot,
#'               'test_identical' returns the test and only the pairs that comprise the same tissues,
#'               'DFtest_identical' returns all the pairs for the test and the identical tissues (as two different objects).
#' @param publish boolean. Default: TRUE. Whether to apply theme_bw()
#'
#' @return the result of the Welch t-test and a figure
#' @export
#'
t_test.CorrMat_identical.vs.different<-function(CorrMat,label.A,label.B,report=TRUE,index="####",
                                                output='DFtest_identical',publish=TRUE){

  if(report) cat(paste0("\n\n",index," Welch Two sample test between the identical tissues and different tissues of ",label.A,
                        " and ",label.B,"\n\n"))
  DFtest <-Reduce(rbind,lapply(colnames(CorrMat),function(x){
    rbind(cbind(DF1=x,
                DF2=reshape2::melt(CorrMat[grep(x,rownames(CorrMat),ignore.case=TRUE),
                                           colnames(CorrMat)[grep(x,colnames(CorrMat),ignore.case=TRUE)],drop=FALSE]),
                type='Same-tissue pairs'),
          cbind(DF1=x,
                DF2=reshape2::melt(CorrMat[grep(x,rownames(CorrMat),ignore.case=TRUE),
                                           !colnames(CorrMat) %in% colnames(CorrMat)[grep(x,colnames(CorrMat),ignore.case=TRUE)]]),
                type='Different tissues pairs')
    )
  }))
  colnames(DFtest)<-c('DF1','DF2','Correlation','Type')
  identical<-DFtest[DFtest$Type=="Same-tissue pairs",]
  different<-DFtest[DFtest$Type=="Different tissues pairs",]
  test<-t.test(identical$Correlation,different$Correlation)

  identical<-identical[order(identical$Correlation,decreasing=TRUE),]
  estimateDF=data.frame(Type=c('Same-tissue pairs','Different tissues pairs'),
                        Correlation=test$estimate)
  if(report) cat.html(test)

  p<-ggplot2::ggplot(DFtest,aes(Type,Correlation))+geom_boxplot()+ annotate('text',label=paste('p-value=',
                                                                                               signif(test$p.value,digits=8)),
                                                                            x=0.75,y=0.1)+geom_text(data=estimateDF,
                                                                                                    aes(label=paste0('m= ',signif(Correlation,digits=4)),
                                                                                                        y=0.02),size=10)# +ylab(paste(simpleCap(METHODCOR),'Correlation'))
  if(publish) p<- p +theme_bw()

  p <- p + ggplot2::theme(axis.text=element_text(size=15),
                          axis.title=element_text(size=15),
                          legend.position='none')

  print(p)
  switch(output,
         'identical'= return(identical),
         'different'= return(different),
         'DFtest'= return(DFtest),
         'test' = return(test),
         'plot' = return(p),
         'test_identical'= return(list(test,identical)),
         'DFtest_identical'= return(list(DFtest,identical)))

}


#' Apply t.test.CorrMat_identical.vs.different and complete the resulting data.frame
#'
#' @param CorrMat correlation matrix between two studies tissue expression (column correlation)
#' @param label.A character string. Name of the first study on which the correlations are based
#' @param label.B character string. Name of the second study on which the correlations are based
#' @param report boolean. Default: TRUE. Whether to output a text report.
#' @param index character string. For the notebook, allows to create new sections.
#' @param publish boolean. Default: TRUE. Whether to apply theme_bw()
#' @param QuantMet character string. Allows to fill in the quantification of the studies.
#' @param datasets Character string. Allows to fill in the names of the studies.
#' @param CorrelationMethod character string. Fill in the name of the correlation method
#'
#' @return a data.frame
#' @export
#'
completeDF.t_test.CorrMat_identical.vs.different<-function(CorrMat,label.A,label.B,report=TRUE,index="####",
                                                  publish=TRUE,QuantMet,datasets,CorrelationMethod){


  newDF<-t_test.CorrMat_identical.vs.different(CorrMat,label.A=label.A,label.B=label.B,
                                               report=report,index=index,output='DFtest',
                                               publish=publish)
  newDF$Tissue.Nb<-base::ncol(CorrMat)
  newDF$QuantMet<-QuantMet
  newDF$datasets<-datasets
  newDF$CorrelationMethod<-CorrelationMethod

  return(newDF)

}


#' Test with a Welch t-test between the group of similar tissues versus the group of different tissues.
#' Take the diagonal and test it versus the lower triangle of the matrix.
#'
#' @param Mat correlation matrix between two studies tissue expression (column correlation)
#' @param yLab character string. Title of the y-axis
#' @param label.diag character string. Label for the pairs in the diagonal of the matrix.
#' @param label.low.tri character string. Label for the pairs in the lower.triangle of the matrix.
#' @param xLab character string. Title of the x-axis. Default: 'Type'.
#' @param pval boolean. Default: TRUE. Whether to add the p-value on the plot.
#' @param showMean boolean. Default: TRUE. Whether to add the mean of each group on the plot.
#' @param output character string that allows to pick the output;
#'               'identical' for the same pairs of tissues,
#'               'different' for the different pairs of tissues,
#'               'test' for the test,
#'               'plot' for the plot and
#'               'test_identical' for a list comprising the identical tissues and the test.
#' @param publish boolean. Default: TRUE. Whether to apply theme_bw().
#'
#' @return depends on 'output'; it can be a data.frame, a figure or the result of the test.
#' @export
#'
t_test.Mat_identical.vs.different<-function(Mat,yLab,
                                            label.diag='Same-tissue pairs',
                                            label.low.tri='Different tissues pairs',
                                            xLab='Type',pval=TRUE,showMean=TRUE,
                                            output='test',publish=TRUE){

  identical<-diag(as.matrix(Mat))
  different<-as.matrix(Mat)[lower.tri(as.matrix(Mat), diag = FALSE)]
  test<-t.test(identical,different)

  estimateDF=data.frame(Type=c(label.diag,label.low.tri),
                        Value=test$estimate)

  DFtest<-rbind(data.frame(Value=identical,Type=label.diag),
                data.frame(Value=different,Type=label.low.tri))

  p <- ggplot2::ggplot(DFtest,aes(Type,Value))+ggplot2::geom_boxplot()
  if(pval) p <- p + ggplot2::annotate('text',label=paste('p-value=',
                                                         signif(test$p.value,digits=8)),
                                      x=0.75,y=0.1)
  if(showMean) p<- p+ ggplot2::geom_text(data=estimateDF,
                                         aes(label=paste0('m= ',signif(Value,digits=4)),
                                             y=0.02),size=10)
  if(publish) p<- p +theme_bw()

  p<-p+ggplot2::xlab(xLab)+ggplot2::ylab(yLab)
  p<-p+ggplot2::theme(axis.text=element_text(size=15),
                      axis.title=element_text(size=15),
                      legend.position='none')

  print(p)
  switch(output,
         'identical'= return(identical),
         'different'= return(different),
         'test' = return(test),
         'plot' = return(p),
         'test_identical'= return(list(test,identical)))
}


## Analyses based on tissue distance ----------------------

#' Compute Jaccard index
#'
#' @param list1 vector of gene names from a first element to compare
#' @param list2 vector of genes names from a second element to compare
#' @param universeSize positive integer. Number of elements in the universe
#' @param pvalue boolean. Default: FALSE. Whether a p-value should be computed
#' @param indexSignif positive integer. Default: 2. Number of significant digits to add for the Jaccard index on the figure.
#' @param pvalSig positive integer. Default: 3. Number of significant digits to add for the p-value on the figure.
#' @param figure boolean. default: FALSE. Whether a figure showing the overlap should be displayed
#' @param Col Colour palette: vector of two character strings. Colours to use to draw the two elements being compared in the venn diagram
#' @param categories vector of two character strings that are the names of the two elements being compared.
#' @param condition character string. Name of the condition compared between the two elements.
#' @param policeFont character string. Name of the font to use.
#' @param saveFile path (character string). Default: NA. Allows to save the figure (if figure=TRUE) at the given path.
#' @param outputType character string. Whether "index" to return the Jaccard index;
#'                                     "pval" for the p-value associated or
#'                                     "list" for a list comprising the index and the pvalue,
#'                                     "intersect" intersection of the two lists.
#' @param ... other parameters for grDevices::cairo_pdf
#'
#'
#' @return depends on outputType. Refer to that argument for the possible output.
#' @export
#'

jaccardInd<-function(list1,list2,universeSize,pvalue=FALSE,
                     indexSignif=2, pvalSig=3,
                     figure=FALSE,Col=c('coral4','darkorchid1'),
                     categories,condition,policeFont='Linux Libertine',
                     saveFile=NA,outputType="index",
                     ...){

  Inter<-length(intersect(list1,list2))
  Uni<-length(union(list1,list2))
  NbPossibleSuccesses<-length(list1)
  Nbdraws<-length(list2)

  jInd<-Inter/Uni

  if(figure)  pvalText<-paste('\nJaccard index = ',signif(jInd,digits=indexSignif))

  if (outputType %in% c("pval","list")) pvalue<-TRUE

  if(pvalue&&missing(universeSize)) {
    warning("Impossible to compute a p-value as the universe size is missing")
    pvalue<-FALSE
    outputType<-"index"
  }

  if(pvalue){
    #p-value to observe Inter or more (thus Inter-1)
    pval<-phyper(Inter-1,NbPossibleSuccesses,universeSize-NbPossibleSuccesses,Nbdraws,lower.tail = FALSE)
    if(figure)  pvalText<-paste(pvalText,'\np-value =',signif(pval,digits=pvalSig))
  }

  if(figure){
    p<-VennDiagram::draw.pairwise.venn(ind=FALSE,
                          fontfamily=policeFont,
                          area1=NbPossibleSuccesses,
                          area2=Nbdraws,
                          cross.area=Inter,
                          category = c(categories[1],categories[2]),
                          cat.col= c(Col[1],Col[2]),
                          cat.pos = c(340,27),
                          cat.cex=rep(2,2),
                          cat.dist=c(0.04,0.04),
                          col = c(Col[1],Col[2]),
                          fill = c(Col[1],Col[2]),
                          cex=rep(2,3),
                          alpha = rep(0.4, 2),
                          euler.d=TRUE,
                          margin=0,
                          scaled=TRUE
    )
    grid.arrange(grobs=list(title=textGrob(condition,gp=gpar(fontsize=50)),
                            textGrob(pvalText,gp=gpar(fontsize=15)),gTree(children=p)),heights=c(1,0.5,5))

    if(!is.na(saveFile)){
      cairo_pdf(filename = saveFile, family= policeFont,...)
      grid.arrange(grobs=list(title=textGrob(condition,gp=gpar(fontsize=50)),
                              textGrob(pvalText,gp=gpar(fontsize=15)),gTree(children=p)),heights=c(1,0.5,5))
      dev.off()
    }

  }

  switch(outputType,
         "index"    = return(jInd),
         "pval"     = return(pval),
         "list"     = return(list("index"=jInd,"pval"=pval)),
         "intersect"= return(intersect(list1,list2))
         )
}

#' Overlap of the tissue-specific proteins with the most specific transcripts to that tissue.
#'
#' @param DFprot numeric data.frame of the protein expression data
#' @param DFmRNA numeric data.frame of the transcriptomic expression data
#' @param pretreatment boolean. Default: FALSE. Whether the considered tissues (columns) and genes (rows)
#'                     need to be filtered to the common set.
#'                     The function requires the same set of tissues/genes to be meaningful.
#'                     This step can be bypassed if it already happened previously.
#' @param mRNAmethod a function name. The function should allow to order the mRNAs based on their specificity rank.
#' @param cutoffProt numeric. Threshold to consider a protein is expressed.
#' @param complete.matrix boolean. Default:TRUE.
#'                        Whether the comparison between the two studies should be done for their whole data.frame
#' @param pvalue boolean. Default: TRUE.
#'               Whether the p-value of the jaccard index should be computed
#' @param ratio numeric (>0). Multiplying factor that
#'              allows increasing the number of transcripts compared with the proteins.
#' @param xprot character string. Default: NA. Needed when complete.matrix is false.
#'              Name of the tissue to consider for the proteomics.
#' @param yrna character string. Default: NA. Needed when complete.matrix is false.
#'              Name of the tissue to consider for the transcriptomics.
#' @param figure boolean. Default: TRUE. Whether the venn diagram should be drawn
#' @param ... more arguments for grDevices::cairo_pdf
#'
#' @return overlap between the proteomics and transcriptomics
#' @export
#'
overlapSpeProtmRNA<-function(DFprot,DFmRNA,pretreatment=FALSE,
                             mRNAmethod='ratioSpe',
                             cutoffProt=0,
                             complete.matrix=TRUE,
                             pvalue=TRUE,
                             ratio=1,
                             xprot=NA,
                             yrna=NA,
                             figure=TRUE,...){
  if(pretreatment){
    commonCond<-intersect(colnames(DFprot),colnames(DFmRNA))
    DFprot<-strip(DFprot[,commonCond])
    DFmRNA<-strip(DFmRNA[,commonCond])

    commonRows<-intersect(rownames(DFprot),rownames(DFmRNA))
    DFprot<-DFprot[commonRows,]
    DFmRNA<-DFmRNA[commonRows,]
  }

  nbGenes<-nrow(DFmRNA)
  #extract the uniquely expressed proteins
  DFprot<-computeBreadth(DFprot,omit.zero=TRUE,
                         threshold=cutoffProt,
                         typeR='unique_tissueList')

  #compute the specificity of RNA in each tissue
  DFmRNA<-eval(call(name=mRNAmethod,DFmRNA))
  conditionCnter<-colnames(DFmRNA)

  if(complete.matrix){
    if(!pvalue){
      resDF<-Reduce(rbind,
                    lapply(setNames(conditionCnter,conditionCnter),function(x){
                      #print(paste("x",x))
                      uniqueProt<-rownames(DFprot[DFprot$condition==x,])
                      #print("uniqueProt")
                      #print(dput(uniqueProt))
                      lg<-length(uniqueProt)
                      limit<-round(ratio*lg)
                      tmp<-as.data.frame(lapply(setNames(conditionCnter,conditionCnter), function(y){
                        #print(paste("y",y))
                        list2comp<-rownames(DFmRNA[order(DFmRNA[,y],decreasing = TRUE),])[1:limit]
                        #print("list2comp")
                        #print(dput(list2comp))
                        jaccardInd(uniqueProt,list2comp,universeSize=nbGenes,pvalue=FALSE,figure=FALSE)
                      }))
                    }))
      rownames(resDF)<-conditionCnter
    }else{
      resDF<-Reduce(rbind,lapply(setNames(conditionCnter,conditionCnter),
                                 function(x){
                                   uniqueProt<-rownames(DFprot[DFprot$condition==x,])
                                   lg<-length(uniqueProt)
                                   limit<-ratio*lg
                                   tmp<-as.data.frame(lapply(setNames(conditionCnter,conditionCnter), function(y){
                                     list2comp<-rownames(DFmRNA[order(DFmRNA[,y],decreasing = TRUE),])[1:limit]
                                     jaccardInd(uniqueProt,list2comp,universeSize=nbGenes,pvalue=TRUE,figure=FALSE,
                                                outputType="pval")
                                   }))
                                 }))
      rownames(resDF)<-conditionCnter
    }
    #row results correspond to the protein tissue samples
    #column results correspond to the mRNA tissue samples
    #example: if column[3]=='Heart' and row[4]=='Liver
    #the comparison is between the transcriptome of the Heart and the Liver proteome.
    return(resDF)
  }else{
    if(missing(yrna)) yrna<-xprot
    uniqueProt<-rownames(DFprot[DFprot$condition==xprot,])
    lg<-length(uniqueProt)
    limit<-ratio*lg
    list2comp<-rownames(DFmRNA[order(DFmRNA[,yrna],decreasing = TRUE),])[1:limit]
    if(figure){
      if(xprot!=yrna) {
        condition=paste(xprot,'~',yrna)
      }else{
        condition=xprot
      }
    }
    jaccardInd(uniqueProt,list2comp,universeSize=nbGenes,pvalue=pvalue,
               figure=figure,condition=condition,...)
  }
}



#' Compute the distance between the tissues within a study.
#' @description the distance is based on the number of genes that are uniquely shared by a number of given tissues only.
#'
#' @param DF numeric data.frame that contains
#' @param nbOfTissues integer (>0) number of tissues in which the genes have to be expressed to be considered
#' @param ColN vector of character strings. Names of the tissues to include in the analysis.
#' @param gid boolean. Default: FALSE. Whether to output the name of the concerned genes instead of their number.
#' @param plot boolean. Default: FALSE. Whether the heatmap should be printed or not
#' @param method character string, name of the agglomeration method to be used (hclust)
#'               This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'               "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#'               "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param notecol character string. Name of the colour to use for the annotation within the heatmap.
#' @param col function to be used to define the colour palette for (gplots) heatmap.2
#'            default: colorRampPalette(c('aliceblue','darkcyan'))
#' @param threshold numeric. Default: 0
#' @param ... for labelling purposes
#'
#' @return a data.frame
#' @export
#'
tissuesDistance<-function(DF,nbOfTissues=2,ColN,gid=FALSE,plot=FALSE,
                          method='ward.D',notecol='gray37',
                          col=colorRampPalette(c('aliceblue','darkcyan')),
                          threshold=0,
                          ...){

  if(missing(ColN)) {
    if(colnames(DF)[ncol(DF)]=='nb.tissues'){
      ColN<-colnames(DF)[-ncol(DF)]
      DF.cleaned<-DF[DF$nb.tissues ==nbOfTissues,]
      DF.cleaned<-DF.cleaned[,ColN]
    }else{
      ColN<-colnames(DF)
      DF.cleaned<-DF[rowSums(DF>threshold)==nbOfTissues,]
    }}

  TissuesCombo<-data.frame(utils::combn(ColN,nbOfTissues,simplify=FALSE),stringsAsFactors=FALSE)

  res<-Reduce(rbind,lapply(TissuesCombo,function(x) {
    x<-as.character(x)
    # we consider rows that were the first element of x is true
    # as such need to make sure that there are indeed rows before the next step:
    if(nrow(DF.cleaned[DF.cleaned[,x[1]] >0,])>0){
      DF.tmp<-DF.cleaned[DF.cleaned[,x[1]] >0,]
      if(!gid){
        count<-nrow(DF.tmp[DF.tmp[,x[2]] >0,])
        #same thing with the second tissues
      }else{
        if(length(rownames(DF.tmp[DF.tmp[,x[2]] >0,]))>0){
          id<-Reduce(function(...) paste(...,sep=";"),rownames(DF.tmp[DF.tmp[,x[2]] >0,]))
        }else{
          id<-'none'
        }
      }
    }else{
      if(!gid){
        count<-0
      }else{
        id<-'none'
      }
    }
    if(!gid){
      return(data.frame(T1=x[1],T2=x[2],count=count))
    }else{
      return(data.frame(T1=x[1],T2=x[2],geneID=id))
    }
  })
  )

  SameTissue<-sapply(colnames(DF.cleaned),function(x){
    if("nb.tissues" %in% colnames(DF)){
      return(setNames(sum(DF[DF$nb.tissues==1,x]),x))
    }else{
      uniquelyExpressed<-DF[rowSums(DF>threshold)==1,]
      return(setNames(sum(uniquelyExpressed[,x]>threshold),x))
    }
  })

  if(!gid){
    res2<-data.frame(T1=res[,2],T2=res[,1],count=res[,3])
    identity<-data.frame(T1=colnames(DF.cleaned),T2=colnames(DF.cleaned),count=SameTissue)
    RES<-rbind(res,identity,res2)

    RES<-dcast(RES,T1 ~ T2,value.var='count')
    rownames(RES)<-RES$T1
    RES<-RES[,-1]

    RES<-RES[sort(rownames(RES)),]
    RES<-RES[,sort(colnames(RES))]

  }else{
    RES<-data.frame(T1=res[,2],T2=res[,1],geneID=res[,3])
  }

  ResToplot<-as.matrix(1/RES)
  diag(ResToplot)<-0
  ResToplot[ResToplot=='Inf']<-1

  if(plot){
    p<-heatmap.2(1-ResToplot,
                 hclustfun=function(x) hclust(x,method = method),
                 trace="none",margins = c(10,10),
                 key=FALSE,
                 #scale='both',
                 distfun=function(c) as.dist(1 - c),
                 cellnote=RES, notecol=notecol,
                 col=col,cexRow =1.4,cexCol = 1.4,
                 revC=TRUE,srtCol=45)
    p
  }

  return(RES)
}


#' Plot heatmap based on the result of tissueDistance
#'
#' @param DF numeric data.frame; output of tissueDistance
#' @param method character string, name of the agglomeration method to be used (hclust)
#'               This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'               "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#'               "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param notecol character string. Name of the colour to use for the annotation within the heatmap.
#' @param col function to be used to define the colour palette for (gplots) heatmap.2
#'            default: colorRampPalette(c('aliceblue','darkcyan'))
#'
#' @return a heatmap (based on (gplots) heatmap.2)
#' @export
#'
plotTissuesDistance<-function(DF,method='ward.D',notecol='gray37',
                              col=colorRampPalette(c('aliceblue','darkcyan'))){

  Toplot<-as.matrix(1/DF)
  diag(Toplot)<-0
  Toplot[Toplot=='Inf']<-1

  heatmap.2(1-Toplot,
            hclustfun=function(x) hclust(x,method = method),
            trace="none",margins = c(10,10),
            key=FALSE,
            #scale='both',
            distfun=function(c) as.dist(1 - c),
            cellnote=as.matrix(DF), notecol=notecol,
            col=col,cexRow =1.4,cexCol = 1.4,
            revC=TRUE,srtCol=45)

}


#' Extract the dendrogram from the heatmap based on tissueDistance
#'
#' @param DF numeric data.frame; output of tissueDistance
#' @param method character string, name of the agglomeration method to be used (hclust)
#'               This should be (an unambiguous abbreviation of) one of "ward.D", "ward.D2", "single",
#'               "complete", "average" (= UPGMA), "mcquitty" (= WPGMA),
#'               "median" (= WPGMC) or "centroid" (= UPGMC).
#' @param row boolean. Default: TRUE. Whether to extract the dendrogram associated to the rows of the heatmap.
#'
#' @return a dendrogram
#' @export
#'
extractDendroFromPlotTD<-function(DF,method='ward.D',row=TRUE){

  Toplot<-as.matrix(1/DF)
  diag(Toplot)<-0
  Toplot[Toplot=='Inf']<-1

  tmp<-heatmap.2(1-Toplot,
                 hclustfun=function(x) hclust(x,method = method),
                 trace="none",
                 key=FALSE,
                 distfun=function(c) as.dist(1 - c),
                 revC=TRUE)

  if(row){
    return(tmp$rowDendrogram)
  }else{
    return(tmp$colDendrogram)
  }
}



#' Create a dendrogram from a tissueDistance data.frame
#'
#' @param DF1 numeric data.frame; output of tissueDistance
#' @param trans boolean. Default: TRUE
#' @param plot boolean. Default: TRUE
#' @param type a character string specifying the type of phylogeny to be drawn;
#'             Default: "radial"
#'             can also be "phylogram", "cladogram", "fan", "unrooted".
#'
#' @return a dendrogram (phylo object)
#' @export
#'
extractDendroTissueDistance<-function(DF1,trans=TRUE,plot=TRUE,type){
  if(trans){
    #dend.DF1<-nj(as.dist(1-DF1))
    dend.DF1<-ape::as.phylo(hclust(as.dist(revCustomMatrixDist(DF1)),method='ward.D'))
  }
  if(plot){
    if(missing(type)){
      print("Default type value: radial.")
      print("Other options: phylogram, cladogram, fan, unrooted")
      type<-'radial'
    }
    plot(dend.DF1,type=type)
  }
  return(dend.DF1)
}


# Misc and depricated ----------------------------------------------------------
# might contain problems in the designed purposes or anything.


#' Calculate the ratio of standard deviation between two data.frames
#'
#' @param DF1 numeric data.frame
#' @param DF2 numeric data.frame
#' @param option1 boolean. Default TRUE. Whether to use the first option to compute the ratio sd:
#' first option: log2(DF1)-log2(DF2); other option: log2(DF1)/log2(DF2)
#' @param pseudocount numeric. Default: 1. What to add to avoid problems with log2(0)
#'
#'
#' @return a numeric data.frame
#' @export
#'
comp_ratio_sd<-function(DF1,DF2,option1=TRUE, pseudocount=1){
  if(option1){
    DF<-log2(DF1+pseudocount)-log2(DF2+pseudocount)
  }else{
    DF<-log2(DF1+pseudocount)/log2(DF2+pseudocount)
  }
  DF$average<-rowMeans(DF)
  DF$sd<-apply(DF[,-ncol(DF)],1,sd)
  DF$avg.sd<-DF$average/DF$sd

  DF<-DF[order(DF$avg.sd,decreasing=TRUE),]

  return(DF)
}


#' Compute the correlation between two studies after their transposition
#' and their transformation with log2(x+1)
#'
#' @param df1 numeric data.frame for the first study
#' @param df2 numeric data.frame for the second study
#'
#' @return a named vector with the name of the first study
#'         and the correlation between the two datasets.
#' @export
#'
comp_cor_log2<-function(df1,df2){
  tdf1<-t(df1)
  tdf2<-t(df2)
  sapply(colnames(tdf1),function(x){
    cor(log2(tdf1[,x] +1),log2(tdf2[,x] +1))
  })
}

# (Depricated) Analyses based on the top genes -------------------

#' Gives the number or the refenrences of shared genes
#' by both studies out of a given number of highest genes to query.
#'
#' @param a data.frame with expression data for the first study
#' @param b data.frame with expression data for the second study
#' @param top positive integer. Number of genes to consider for the analysis;
#'            if 0 (default), all genes are considered.
#' @param filta character vector or positive integer vector.
#'              rownames or indices of the rows to consider for the first study.
#' @param filtb character vector or positive integer vector.
#'              rownames or indices of the rows to consider for the second study.
#' @param cutoff numeric. cutoff under which the expression is considered as noise
#' @param LG Boolean, topCommonGenes return the length if TRUE, else the list of genes is returned
#'
#' @return depends on LG. Number or the list of genes that are shared for the highest expressed
#' @export
#'
topCommonGenes <- function(a,b,top=0,filta,filtb,cutoff,LG=TRUE){

  if(!missing(cutoff)){#only the genes expressed at least once above cutoff are kept
    a<-df.cutoff(a,cutoff)
    b<-df.cutoff(b,cutoff)
  }

  if (top==0){#whole datasets are considered
    common.genes<-base::intersect(row.names(a),row.names(b))
  }
  else# only the 'top' genes are considered
  {
    a.order<-a[base::order(filta,decreasing=TRUE),]
    b.order<-b[base::order(filtb,decreasing=TRUE),]
    if (length(row.names(a.order))>=top){
      a.ind<-row.names(a.order)[1:top]
    }else{
      a.ind<-row.names(a.order)
      message("The first dataset has less genes than the given number of top genes to be considered")
    }
    if (length(row.names(b.order))>=top){
      b.ind<-row.names(b.order)[1:top]
    }else{
      b.ind<-row.names(b.order)
      message("The second dataset has less genes than the given number of top genes to be considered")
    }
    common.genes<-base::intersect(a.ind,b.ind)
  }
  if (LG==TRUE) {
    return(length(common.genes))
  }else{
    return(common.genes)
  }
}

#' Compute the correlation matrix for the correlation of the two studies' tissues (columns)
#'
#' @param a data.frame with expression data for the first study
#' @param b data.frame with expression data for the second study
#' @param filta character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the first study.
#' @param same boolean. Default:FALSE. Whether the tissues (columns)
#'             from the second study should be identical to the first study.
#' @param filtb character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the second study.
#' @param top  positive integer. Number of genes to consider for the analysis;
#'            if 0 (default), all genes are considered.
#' @param use character string expliciting how missing values should be handled.
#'            can be either  "everything", "all.obs", "complete.obs", "na.or.complete",
#'            or "pairwise.complete.obs" (default).
#'            To fasten the process if there is no missing value, use "everything".
#' @param method character string. Method for base::cor, either "pearson", "spearman" or
#'              "kendall".
#' @param nameA character string to use for identifying the first study
#' @param nameB character string to use for identifying the second study
#'
#' @return a matrix with the correlation of the columns of the two data.frames
#' @export
#'
topgenes.corCluster <- function(a,b,filta=NA,same=FALSE,filtb=NA,top=0,
                                use="pairwise.complete.obs",method="pearson",
                                nameA=deparse(substitute(a)),nameB=deparse(substitute(b)))
{


  if (top==0){
    common.genes<-base::intersect(row.names(a),row.names(b))
  }else{

    if (!is.na(filta[1])){	## <==> if a filter has been provided
      a.order<-a[base::order(a[,filta],decreasing=TRUE),]
    }else{
      a.tmp<-a[base::order(a,decreasing=TRUE),]
      a.order<-na.omit(a.tmp)#some lines with NA are added, they're deleted
    }
    if(same){
      filtb<-filta
      printDebug(filtb)
    }
    if (!is.na(filtb[1])){
      ## <==> filtb has been provided (call or with 'same')
      b.order<-b[base::order(b[,filtb],decreasing=TRUE),]
    }else{
      b.tmp<-b[base::order(b,decreasing=TRUE),]
      b.order<-na.omit(b.tmp)#some lines with NA are added, they're deleted
    }
    a.ind <- row.names(a.order)[1:top]
    b.ind <- row.names(b.order)[1:top]
    common.genes<-base::intersect(a.ind,b.ind)
  }


  if (length(common.genes)>2) {

    #creation of one dataset (with all the data)
    new.all<-cross.selectCond(a,b,name1=nameA,name2=nameB)
    #calcul of the correlation
    data.cor=cor(new.all,use=use, method=method)
    return(data.cor)

  }else #there is no genes in common
  {
    vide<-matrix(data=NA,nrow=length(colnames(a)),ncol=length(colnames(b)),dimnames=list(colnames(a),colnames(b)))
    return(vide)
  }
}


#' Draw (after computation) the correlation matrix for the correlation of the two studies' tissues (columns)
#'
#' @param a data.frame with expression data for the first study
#' @param b data.frame with expression data for the second study
#' @param filta character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the first study.
#' @param same boolean. Default:FALSE. Whether the tissues (columns)
#'             from the second study should be identical to the first study.
#' @param filtb character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the second study.
#' @param top positive integer. Number of genes to consider for the analysis;
#'            if 0 (default), all genes are considered.
#' @param use character string expliciting how missing values should be handled.
#'            can be either  "everything", "all.obs", "complete.obs", "na.or.complete",
#'            or "pairwise.complete.obs" (default).
#'            To fasten the process if there is no missing value, use "everything".
#' @param method character string. Method for base::cor, either "pearson", "spearman" or
#'              "kendall".
#' @param nameA character string to use for identifying the first study
#' @param nameB character string to use for identifying the second study
#' @param png boolean. Default:TRUE. Whether to save the figure as a png
#'            instead of the standard output
#'
#' @return a figure, either on the standard output or as png
#' @export
#'
topgenes.corClusterDraw <- function(a,b,filta=NA,same=FALSE,filtb=NA,
                                    top=0,use="pairwise.complete.obs",
                                    method="pearson",
                                    nameA=deparse(substitute(a)),
                                    nameB=deparse(substitute(b)),png=TRUE){

  new=topgenes.corCluster(a,b,filta=filta,same=same,filtb=filtb,top=top,nameA=nameA,nameB=nameB,method=method,use=use)
  col=colorRampPalette(c('aliceblue','darkcyan'))
  if (png){
    tmp=paste(nameA,'_',sep="") #initialization of the filename
    if(!is.na(filta)){
      tmp=paste(tmp,filta,sep="")
    }
    tmp=paste(tmp,paste('_',nameB,sep=""),sep="")
    if (same) filtb <- filta
    if(!is.na(filtb)){
      tmp<-paste(tmp,paste('_',filtb,sep=""),sep="")
    }
    tmp<-paste(tmp,'_',sep="")
    png(file=paste(tmp,paste(method,"_CorClust.png",sep=""),sep=""), bg="transparent")
    graphics::par(oma=c(2.5,2.5,2.5,2.5))
    heatmap.2(new, distfun=function(c) as.dist(1 - c), trace="none", dendrogram="row", col=col)
    dev.off()
  }else{
    graphics::par(oma=c(2.5,2.5,2.5,2.5))
    heatmap.2(new, distfun=function(c) as.dist(1 - c), trace="none", dendrogram="row", col=col)
  }
  return(new)
}


#' Compute the correlation matrix for the correlation of the two studies' tissues (columns)
#' after they have been transformed with log
#'
#' @param a data.frame with expression data for the first study
#' @param b data.frame with expression data for the second study
#' @param filta character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the first study.
#' @param same boolean. Default:FALSE. Whether the tissues (columns)
#'             from the second study should be identical to the first study.
#' @param filtb character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the second study.
#' @param top positive integer. Number of genes to consider for the analysis;
#'            if 0 (default), all genes are considered.
#' @param use character string expliciting how missing values should be handled.
#'            can be either  "everything", "all.obs", "complete.obs", "na.or.complete",
#'            or "pairwise.complete.obs" (default).
#'            To fasten the process if there is no missing value, use "everything".
#' @param method character string. Method for base::cor, either "pearson", "spearman" or
#'              "kendall".
#' @param nameA character string to use for identifying the first study
#' @param nameB character string to use for identifying the second study
#' @param cleanInfinite boolean. Default: TRUE. whether to use clean.infinite to handle log(0)
#' @param pseudocount numeric. Default: 1. Pseudocount in case for log(x+pseudocount) when clean.infinite is not used for log(0)
#' @return a matrix with the correlation of the columns of the two data.frames once log2-transformed
#' @export
#'
topgenes.corCluster_log<- function(a,b,filta=NA,same=FALSE,filtb=NA,top=0,
                                   use="pairwise.complete.obs",method="pearson",
                                   nameA=deparse(substitute(a)),
                                   nameB=deparse(substitute(b)),
                                   cleanInfinite=TRUE,
                                   pseudocount=1){

  if (top==0){
    common.genes<-base::intersect(row.names(a),row.names(b))
  }else{

    if (!is.na(filta[1])){	## <==> if a filter has been provided
      a.order<-a[base::order(a[,filta],decreasing=TRUE),]
    }else{
      a.tmp<-a[base::order(a,decreasing=TRUE),]
      a.order<-na.omit(a.tmp)#some lines with NA are added, they're deleted
    }
    if(same){
      filtb<-filta
      printDebug(filtb)
    }
    if (!is.na(filtb[1])){
      ## <==> filtb has been provided (call or with 'same')
      b.order<-b[base::order(b[,filtb],decreasing=TRUE),]
    }else{
      b.tmp<-b[base::order(b,decreasing=TRUE),]
      b.order<-na.omit(b.tmp)#some lines with NA are added, they're deleted
    }
    a.ind<-row.names(a.order)[1:top]
    b.ind<-row.names(b.order)[1:top]
    common.genes<-base::intersect(a.ind,b.ind)
  }


  if (length(common.genes)>2) {

    #creation of one dataset (with all the data)
    new.all<-cross.selectCond(a,b,name1=nameA,name2=nameB)
    #calcul of the correlation
    if(cleanInfinite){
      data.cor<-cor(clean.infinite(log(new.all)),use=use, method=method)
    }else{
      data.cor<-cor(log(new.all+pseudocount),use=use, method=method)
    }
    return(data.cor)

  }else #there is no genes in common
  {
    vide<-matrix(data=NA,nrow=length(colnames(a)),ncol=length(colnames(b)),dimnames=list(colnames(a),colnames(b)))
    return(vide)
  }
}


#' Draw (after computation) the correlation matrix for the correlation of the two studies' tissues (columns)
#' after they have been transformed with log
#'
#' @param a data.frame with expression data for the first study
#' @param b data.frame with expression data for the second study
#' @param filta character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the first study.
#' @param same boolean. Default:FALSE. Whether the tissues (columns)
#'             from the second study should be identical to the first study.
#' @param filtb character vector or positive integer vector.
#'              colnames or indices of the columns to consider for the second study.
#' @param top  positive integer. Number of genes to consider for the analysis;
#'            if 0 (default), all genes are considered.
#' @param use character string expliciting how missing values should be handled.
#'            can be either  "everything", "all.obs", "complete.obs", "na.or.complete",
#'            or "pairwise.complete.obs" (default).
#'            To fasten the process if there is no missing value, use "everything".
#' @param method character string. Method for base::cor, either "pearson", "spearman" or
#'              "kendall".
#' @param nameA character string to use for identifying the first study
#' @param nameB character string to use for identifying the second study
#' @param cleanInfinite boolean. Default: TRUE. whether to use clean.infinite to handle log(0)
#' @param pseudocount numeric. Default: 1. Pseudocount in case for log(x+pseudocount) when clean.infinite is not used for log(0)
#' @param png boolean. Default:TRUE. Whether to save the figure as a png
#'            instead of the standard output
#'
#' @return a heatmap of the correlation matrix on the log-transformed expression data of two studies
#' @export
#'
topgenes.corClusterDraw_log<-function(a,b,filta=NA,same=FALSE,filtb=NA,top=0,
                                      use="pairwise.complete.obs",method="pearson",
                                      nameA=deparse(substitute(a)),nameB=deparse(substitute(b)),
                                      png=TRUE,
                                      cleanInfinite=TRUE,
                                      pseudocount=1){

  new=topgenes.corCluster_log(a,b,filta=filta,same=same,filtb=filtb,top=top,nameA=nameA,nameB=nameB,method=method,use=use)
  col=colorRampPalette(c('aliceblue','darkcyan'))
  if (png){
    tmp=paste(nameA,'_',sep="") #initialization of the filename
    if(!is.na(filta)){
      tmp=paste(tmp,filta,sep="")
    }#endif !is.na(filta)
    tmp=paste(tmp,paste('_',nameB,sep=""),sep="")
    if (same) filtb <- filta
    if(!is.na(filtb)){
      tmp<-paste(tmp,paste('_',filtb,sep=""),sep="")
    }
    tmp<-paste(tmp,'_',sep="")
    png(file=paste(tmp,paste(method,"_CorClust_log2.png",sep=""),sep=""), bg="transparent")
    graphics::par(oma=c(2.5,2.5,2.5,2.5))
    heatmap.2(new, distfun=function(c) as.dist(1 - c), trace="none", dendrogram="row", col=col)
    dev.off()
  }else{
    graphics::par(oma=c(2.5,2.5,2.5,2.5))
    heatmap.2(new, distfun=function(c) as.dist(1 - c), trace="none", dendrogram="row", col=col)
  }
  return(new)
}


# Based on the breadth of expression of the genes as a first step

#' Create matrix of jaccard indices based on the output of ExpressedGeneInTissue
#'
#' @param DFa data.frame
#' @param geneList list of the genes found in tissue
#' @param ... other parameters for ExpressedGeneInTissue
#'
#' @return a matrix
#' @export
#'
internalJaccardIndMatrix<-function(DFa,geneList, ...){
  if(missing(geneList)) geneList<-ExpressedGeneInTissue(DFa, ...)
  res1<-as.data.frame(t(utils::combn(colnames(DFa),2)),stringsAsFactors = FALSE)
  res1$val<-sapply(1:nrow(res1),function(x){
    return(jaccardInd(unlist(geneList[res1[x,1]]),unlist(geneList[res1[x,2]])))
    })
  res2<-data.frame("V1"=res1$V2,"V2"=res1$V1,val=res1$val,stringsAsFactors = FALSE)
  res3<-data.frame("V1"=colnames(DFa),"V2"=colnames(DFa),val=1,stringsAsFactors = FALSE)
  res<-rbind(res1,res2,res3)
  return(as.matrix(suppressMessages(reshape2::acast(res, V1~V2))))
}

