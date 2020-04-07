# Correlation based ------------------------------------------------

#' Compute for the genes the correlation of their expression between two datasets
#' across different tissues
#' @description Both data.frames have to share the same number of dimensions.
#'              For meaningful correlation,
#'              the rows and the columns should be identical between the two datasets.
#'
#' @param DF1 numeric data.frame
#' @param DF2 numeric data.frame.
#' @param method a character string indicating which correlation coefficient is to be computed.
#' One of "pearson", "spearman", "kendall" or 'pearson_log2'
#' (for 'pearson_log2' the data is first log2 transformed before the computation of the correlation values)
#'
#' @return a named vector that contains the correlation value of each gene (used as name)
#' @export
#'
compute_gene_cor<-function(DF1,DF2,method){
  #initialise
  nbGenes=nrow(DF1)
  CorVec<-setNames(rep(0,nbGenes),rownames(DF1))

  #compute
  if(method!='pearson_log2'){
    CorVec<-sapply(names(CorVec),function(x){
      cor(as.numeric(DF1[x,]),as.numeric(DF2[x,]),method=method)
    })
  }else{
    DF1<-log2(DF1+1)
    DF2<-log2(DF2+1)
    CorVec<-sapply(names(CorVec),function(x){
      cor(as.numeric(DF1[x,]),as.numeric(DF2[x,]),method='pearson')
    })
  }

  return(CorVec)
}

#' Allows to save to a file the correlation of the genes expression
#' between two datasets across different tissues
#'
#' @param DF1 numeric data.frame
#' @param DF2 numeric data.frame
#' @param method a character string indicating which correlation coefficient
#' is to be computed. #' One of "pearson", "kendall", "spearman".
#' @param save boolean; default: TRUE. Whether the correlation should be saved in a file
#' @param out.df1 character string to use for DF1 in the saved file label
#' @param out.df2 character string to use for DF2 in the saved file label
#' @param ... other possible parameters for \code{\link[barzinePhdR]{saveToFile}}
#'
#' @return numeric named vector while also saving it to a file
#' @export
#'
compute_gene_cor_save<-function(DF1,DF2,method,out.df1,out.df2,save,...){
  CorVec<-compute_gene_cor(DF1,DF2,method=method)
  if (save) message<-saveToFile(CorVec,
                                filename=paste0(method,
                                                '_Correlation_',out.df1,'_VS_',
                                                out.df2,'_realData'),
                                ...)
  message(message)
  return(CorVec)
}


#' Shuffle randomly one expression dataset before computing the correlation with the same labelled rows
#' from another dataset across the same set of different tissues
#'
#' @param x numeric data.frame for a first expression study
#' @param y numeric data.frame for a second expression study
#' @param method a character string indicating which correlation coefficient is to be computed.
#'               One of "pearson", "spearman", "kendall".
#' @return a random permutation for compute_gene_cor
#' @export
#'
compute_permuted_cor<-function(x,y,method="spearman"){
  TMP<-x[base::sample(nrow(x)),] #create a permuted (scramble the rows)
  result<-sapply(1:nrow(y),function(i){
    cor(as.numeric(TMP[i,]),as.numeric(y[i,]),method=method)
  })
  names(result)<-rownames(y)
  return(result)
}


#' Optimised to work with clustermq (to be used for LSF jobs)
#' and shuffle randomly one expression dataset before computing the correlation with the same labelled rows
#' from another dataset across the same set of different tissues
#'
#' @param rep positive integer. The current iteration for
#' @param x numeric data.frame for a first expression study
#' @param y numeric data.frame for a second expression study
#' @param g.names character vector. Names of the genes.
#' @param method a character string indicating which correlation coefficient is to be computed.
#'               One of "pearson", "spearman", "kendall".
#' @param random positive integer. Random seed to add to "rep" so while the shuffling is random,
#'               it can be reproduced.
#'
#' @return a random shuffling instance for compute_gene_cor
#' @export
#'
compute_permuted_corQ<-function(rep,x,y,g.names,method='spearman',random=145){
  #as the data is given as list
  set.seed(rep+random)
  x<-matrix(unlist(x),nrow=length(unlist(g.names)))
  y<-matrix(unlist(y),nrow=length(unlist(g.names)))
  TMP<-x[base::sample(nrow(x)),]
  result<-sapply(1:nrow(y), function(i){
    cor(as.numeric(TMP[i,]),as.numeric(y[i,]),method=method)
  })
  names(result)<-unlist(g.names)
  return(result)
}


#' Wrapper to run the random shuffling of the gene correlation as several jobs on LSF
#'
#' @param DF1 numeric data.frame for the first study
#' @param DF2 numeric data.frame for the second study
#' @param method a character string indicating which correlation coefficient is to be computed.
#'               One of "pearson", "spearman", "kendall".
#' @param out.df1 character string. Label for the first study. Part of the filename.
#' @param out.df2 character string. Label for the second study. Part of the filename.
#' @param repTimes positive integer. Number of randomisation to create. e.g. 10,000
#' @param seed random seed to maintain repeatability
#' @param nJobs positive integer. Number of simultaneous jobs.
#' @param save boolean. Default: TRUE. The function tries to save in any case,
#'             this argument defines if there is already a file with the same name
#'             if it should be overwrite or not.
#' @param importPath path to ebits --- newest versions are found in \href{https://github.com/mschubert/ebits}{M. Schubert's github profile}
#' @param ... other parameters for ebits::hpc or saveToFile
#'
#' @return Split the several randomisation of the genes correlation as different jobs on LSF
#' @export
#'
wrap.compute_permuted_corQ_save<-function(DF1,DF2,
                                          method,out.df1,out.df2,repTimes,seed,
                                          nJobs,save=TRUE,
                                          importPath='~mitra/homeshare/R321/ebits', ...){

  options(import.path=importPath)
  hpc=modules::import('hpc')

  CorRand.all<-data.frame(hpc$Q(compute_permuted_corQ,rep=1:repTimes,x=list(DF1),y=list(DF2),
                                g.names=list(rownames(DF1)),random=seed,n_jobs=nJobs,...))
  colnames(CorRand.all)<-sapply(1:repTimes,function(x) paste0('v',x))
  saveToFile(CorRand.all,filename=paste0(method,'_Correlation_',
                                                   out.df1,'_VS_',out.df2,'_simulData'),
             overwrite=save,...)
  return(CorRand.all)
}



#' Plot the correlation of different studies comparison that are ranked in decreasing order of correlation level.
#'
#' @param cor1 named numeric vector. Comprises all the correlations between two studies for each genes.
#' @param cor2 named numeric vector. Comprises all the correlations between two other studies for each genes.
#' @param ref named numeric vector. Comprises all the correlations between two studies that can be used as reference.
#' @param simul1 named numeric vector or data.frame. Randomisation of cor1
#' @param simul2 named numeric vector or data.frame. Randomisation of cor2
#' @param seed integer. Default: 1235. Allows to always pick the same simulation when data.frames are supplied.
#' @param aggregateBool boolean. Default:FALSE. Whether all the simulations should be aggreagated in one value for each gene per data.frame
#' @param aggregMethod character string for a funtion. Default:'rowMeans'.
#' @param title character string. Title of the plot.
#' @param xtitle character string. Label of the x-axis.
#' @param ytitle character string. Label of the y-axis.
#' @param unit character string. Unit in which the expression of the first study (Proteomics) is quantified
#' @param base_family character string. Name of the font to use.
#' @param base_size numeric. Size of the font to use (as a base)
#' @param legendText numeric. Size of the legend text.
#' @param ColPalette vector of character string to customise
#'                   the colours of "cor1", "cor2","simul1",...
#'
#' @param geneList vector of character string. Gene identifiers that should be included in the figure.
#'                 default: common ids between cor1, cor2, simul1, simul2 and ref
#' @param out character string. Either "point" or "density"
#'
#' @return a figure
#' @export
#'
sortedGenesCorr<-function(cor1,cor2,ref,simul1,simul2,
                          seed=1235,
                          aggregateBool=FALSE,
                          aggregMethod='rowMeans',
                          title='',
                          xtitle='Rank (sorted by decreasing order of correlation)',
                          ytitle='',
                          unit=' (PPKM)',
                          base_family="Linux Libertine",
                          base_size=11,
                          legendText=0.92,
                          ColPalette,
                          geneList,
                          out='point'){

  set.seed(seed)

  if(missing(ColPalette)){
    if(missing(cor2)){
      ColPalette<-setNames(c('plum','grey68','forestgreen'),
                           c("Protein/mRNA pairs",
                             "Randomised Protein/mRNA pairs",
                             "mRNA/mRNA pairs ('ideal' reference)"))
    }else{
      ColPalette<-setNames(c('plum','grey68','plum1','grey81','lightgreen'),
                           c(paste0('Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             "Reference: Uhl\u00e9n mRNA/GTEx mRNA pairs ('ideal' case)"))
    }
  }

  categoriesNames<-names(ColPalette)

  if(!missing(cor2)){
    commonNames<-Intersect(names(cor1),names(cor2),names(ref))
  }else{
    commonNames<-intersect(names(cor1),names(ref))
  }

  if(!missing(geneList)) commonNames<-geneList

  maxSeq<-length(commonNames)
  index<-1:maxSeq

  CorDF<-data.frame(Ranked.cor=index,
                    Cor=sort(cor1[commonNames],decreasing=TRUE),
                    Comparison=categoriesNames[1])

  if(!missing(cor2)){
    CorDF<-rbind(CorDF,
                 data.frame(Ranked.cor=index,
                            Cor=sort(cor2[commonNames],decreasing=TRUE),
                            Comparison=categoriesNames[3]))
  }

  if(!missing(simul1)){
    if(!aggregateBool) {
      CorDF<-rbind(CorDF,
                   data.frame(Ranked.cor=index,
                              Cor=sort(simul1[commonNames,sample(ncol(simul1),1)],decreasing=TRUE),
                              Comparison=categoriesNames[2],stringsAsFactors =FALSE))
    }else{
      simul1<-simul1[commonNames,]
      sim1<-apply(simul1,2,sort)
      sim1<-eval(call(name=aggregMethod,as.matrix(sim1)))
      sim1<-sort(sim1,decreasing = TRUE)
      CorDF<-rbind(CorDF,
                   data.frame(Ranked.cor=index,
                              Cor=sim1,
                              Comparison=categoriesNames[2],stringsAsFactors =FALSE))
    }
  }

  if(!missing(simul2)){
    if(!aggregateBool){
      CorDF<-rbind(CorDF,
                   data.frame(Ranked.cor=index,
                              Cor=sort(simul2[commonNames,sample(ncol(simul2),1)],decreasing=TRUE),
                              Comparison=categoriesNames[4]))
    }else{
      simul2<-simul2[commonNames,]
      sim2<-apply(simul2,2,sort)
      sim2<-eval(call(name=aggregMethod,as.matrix(sim2)))
      sim2<-sort(sim2,decreasing = TRUE)
      CorDF<-rbind(CorDF,
                   data.frame(Ranked.cor=index,
                              Cor=sim2,
                              Comparison=categoriesNames[2],stringsAsFactors =FALSE))
    }
  }
  if(!missing(ref)){
    CorDF<-rbind(CorDF,
                 data.frame(Ranked.cor=index,
                            Cor=sort(ref[commonNames],decreasing=TRUE),
                            Comparison=categoriesNames[length(categoriesNames)]))
  }

  CorDF<-as.data.frame(CorDF)

  if(out=='point'){
    Inter05=sum(CorDF[as.character(CorDF$Comparison)==categoriesNames[1],'Cor']>0.5)
    Inter0=sum(CorDF[as.character(CorDF$Comparison)==categoriesNames[1],'Cor']>0)

    p<-ggplot(CorDF,aes(x=Ranked.cor,y=Cor,group=Comparison,colour=Comparison))
    p<-p+guides(color=guide_legend(nrow=2,byrow=TRUE))+guides(color=guide_legend(ncol=1))+theme(legend.position='bottom')
    p<-p+scale_colour_manual(values=ColPalette)
    p<-p+labs(title = title,x=xtitle,y=ytitle)
    p<-p+geom_vline(colour='steelblue1',xintercept=Inter05)
    p<-p+geom_vline(colour='slategrey',xintercept=Inter0)
    p<-p+geom_hline(colour='royalblue4',yintercept = 0.8)
    p<-p+geom_hline(colour='steelblue1',yintercept=0.5)
    p<-p+geom_hline(colour='slategrey',yintercept=0)
    p<-p+annotate('text',y=-0.9,x=Inter05/2-2,label=paste0('<--',Inter05,' genes -->'))
    p<-p+annotate('text',y=-0.9,x=(Inter0-Inter05)/2+Inter05,label=paste('<--  ',Inter0-Inter05+1,'genes  -->'))
    p<-p+annotate('text',y=-0.9,x=(length(commonNames)-Inter0)/2+Inter0,label=paste('<-- ',length(commonNames)-Inter0,'genes -->'))

    p<-p+theme_classic(base_family = base_family,base_size = base_size)

    p<-p+geom_rangeframe(color='black')+theme(axis.line = element_blank())

    if(!missing(cor2)){
      legend1<-g_legend(p+geom_point()+theme(legend.text= element_text(size = rel(legendText)),legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(ncol=2)))
    }else{
      legend1<-g_legend(p+geom_point()+theme(legend.text= element_text(size = rel(legendText)),legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(title.position = "left")))
    }
    gp1<-p+geom_point(alpha=0.15)+theme(plot.title=element_blank(),legend.position="none")#+geom_rangeframe(color='black')

    return(suppressWarnings(grid.arrange(gp1,legend1,ncol=1,heights=c(5,1))))
  }else{
    if(out=='density'){
      #if(!missing(cor2)){
      #  legend1<-g_legend(p+geom_density()+theme(legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(ncol=2)))
      #}else{
      #  legend1<-g_legend(p+geom_density()+theme(legend.margin=margin(c(0,0,0,1))))
      #}
      #gp1<-p+geom_density(alpha=0.15)+theme(legend.position="none")

      #return(suppressWarnings(grid.arrange(gp1,legend1,ncol=1,heights=c(5,1))))

      p<-ggplot(CorDF,aes(x=Cor,group=Comparison,colour=Comparison))+geom_density()
      p<-p+guides(color=guide_legend(nrow=2,byrow=TRUE))+guides(color=guide_legend(ncol=1))+theme(legend.position='bottom')
      p<-p+scale_colour_manual(values=ColPalette)
      return(p)

    }else{
      return(CorDF)
    }
  }
}


# Intersection based -----------------------------------------------

#' Give the number of the common genes between five datasets at each rank
#' @description The genes need to be ranked based on
#'              a descriptor prior to the use of this function
#'
#' @param A Named numeric vector for the first dataset
#' @param B Named numeric vector for the second dataset
#' @param C Named numeric vector for the third dataset
#' @param D Named numeric vector for the fourth dataset
#' @param E Named numeric vector for the fifth dataset
#' @param decreasing logical, default TRUE.
#'                   Whether the objects should be sorted prior to the comparison
#'
#' @return data.frame ready to be plotted with \code{\link[ggplot2]{ggplot}}-based function
#' @export
#'
overlapSortCumul5DF<-function(A,B,C,D,E,decreasing=TRUE){
  common<-Intersect(names(A),names(B),names(C),names(D),names(E))
  A<-A[common]
  B<-B[common]
  C<-C[common]
  D<-D[common]
  E<-E[common]
  A<-names(sort(A,decreasing=decreasing))
  B<-names(sort(B,decreasing=decreasing))
  C<-names(sort(C,decreasing=decreasing))
  D<-names(sort(D,decreasing=decreasing))
  E<-names(sort(E,decreasing=decreasing))
  res<-data.frame(seq=1:length(common),
                  overlapCount=sapply(1:length(common),
                                      function(x){
                                        length(Intersect(A[1:x],
                                                         B[1:x],
                                                         C[1:x],
                                                         D[1:x],
                                                         E[1:x]))
                                      })
  )
  res$overlap<-res$overlapCount/res$seq
  return(res)
}


#'  Give the number of the common genes between two datasets at each rank
#'
#' @param A Named numeric vector for the first dataset
#' @param B Named numeric vector for the second dataset
#' @param decreasing logical, default TRUE.
#'                   Whether the objects should be sorted prior to the comparison
#'
#' @return data.frame ready to be plotted with \code{\link[ggplot2]{ggplot}}-based function
#' @export
#'
overlapSortCumul2DF<-function(A,B,decreasing=TRUE){
  common<-Intersect(names(A),names(B))
  A<-A[common]
  B<-B[common]
  A<-names(sort(A,decreasing=decreasing))
  B<-names(sort(B,decreasing=decreasing))
  res<-data.frame(seq=1:length(common),
                  overlapCount=sapply(1:length(common),
                                      function(x){
                                        length(Intersect(A[1:x],
                                                         B[1:x]))
                                      })
  )
  res$overlap<-res$overlapCount/res$seq
  return(res)
}

# Based on ranked genes ----

#' Plot the correlation of the gene expression between the same tissues of two studies
#' as a function of the number of genes included in the correlation computation which
#' varies based on the minimal threshold of expression of the genes.
#'
#' @param DFa data.frame that contains the expression data for one tissue across several studies
#' @param sepSeq positive integer vector. Give the indices of the genes for which the expression
#'               should be plotted once ordered.
#' @param step positive integer. The step to use between the rank to plot.
#' @param method a character string indicating which correlation coefficient is to be computed.
#' One of "pearson", "kendall" or "spearman".
#' @param use character string which allows to pick the way missing data are handled.
#'            one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param filename character string. Path where the figure should be saved.
#' @param fontfamily character string, name of the font to use. Default: "Linux Libertine"
#' @param fontSize positive numeric. Size of the font to use as a base.
#' @param simple boolean. Default: FALSE. Whether to simplify the plot.
#' @param sep character string. Default: '('. Allows to retrieve the name of the study from the colnames
#'            (preformated colnames)
#' @param verbose boolean. Default: FALSE. output additional messages.
#' @param position.legend character string. Where to position the legend in regard of the plot.
#'            Default: "bottom", can be also 'right', 'left' and 'top'.
#' @param legendN positive integer. Allows to apply a number of rows to the legend.
#' @param point boolean. Default: TRUE. Whether to draw the curves as a separate point for each value
#' @param smooth boolean. Default:FALSE. Whether to add a smooth line ("lm" method)
#' @param line boolean. Default:FALSE. Whether to draw the curves as lines.
#' @param noTitleLegend boolean. Default: TRUE. Whether to remove the title legend.
#' @param Colorpalette Named vector of character strings that allows to customise the colours in the plot
#' @param logscaleX boolean. Default: FALSE. Whether to apply a log10 transformation to the x-axis scale.
#' @param ylabTitle character string. Allows to redefine the label of the y-axis.
#' @param shorten boolean. Default: FALSE. Allows to shorten the names of the plotted objects (columns)
#' @param Title character string. Title of the plot
#' @param unit character string. Allows to fill in the unit in which the genes are quantified
#'
#' @return a figure
#' @export
#'
evolCorrCumul<-function(DFa,sepSeq,step=10,method='pearson',use='pairwise.complete.obs',
                        filename,fontfamily='Linux Libertine',fontSize=11,simple=FALSE,sep='(',
                        verbose=FALSE,
                        position.legend='bottom',legendN=2,point=TRUE,smooth=FALSE,line=FALSE,
                        noTitleLegend=TRUE,Colorpalette,logscaleX=FALSE,ylabTitle,
                        shorten=FALSE,Title,unit){

  XLAB='Expression level cut-off'
  autopalette<-FALSE
  if(!missing(ylabTitle)){
    newylab<-TRUE
  }else{
    newylab=FALSE
  }

  if(missing(sepSeq)) sepSeq<-c(seq(1,max(DFa),step))

  if(!simple){
    if(shorten){
      colnames(DFa)<-lapply(colnames(DFa),function(x){
        gsub("[\\( \\)]",'',strsplit(x,' ')[[1]][2])
      })
    }
    DFtofill<-as.data.frame(t(data.frame(combn(colnames(DFa),2,simplify=FALSE))),
                            stringsAsFactors = FALSE)
    colnames(DFtofill)<-c('E1','E2')
  }else{
    ColNamesDF<-sort(colnames(DFa))
    DFtofill<-data.frame(E1=ColNamesDF[c(TRUE,FALSE)], #retrieve odd position
                         E2=ColNamesDF[c(FALSE,TRUE)], #retrieve even position
                         stringsAsFactors = FALSE)
  }
  rownames(DFtofill)<-1:nrow(DFtofill)

  if(length(sepSeq)<2){
    message("sequence of breaks too small a random value has been picked")
    sepSeq[2]<-sample(DFa[,1],1)
  }

  DFtofill2<-DFtofill
  DFtofill2$Cut<-sepSeq[1]
  for(i in sepSeq[2:length(sepSeq)]){
    DFtoFillTemp<-DFtofill
    DFtoFillTemp$Cut<-i
    DFtofill2<-rbind(DFtofill2,
                     DFtoFillTemp)
  }

  DFtofill2$Correlation=sapply(1:nrow(DFtofill2),function(x){
    if(verbose) print(paste('Processing',DFtofill2[x,c('E1','Cut')],'...'))
    temp<-strip(DFa[,
                    unlist(DFtofill2[x,c('E1','E2')])
                    ],
                'ge',DFtofill2[x,'Cut']
    )
    if(nrow(temp)>0){
      return(cor(temp[,1],temp[,2],method=method,use=use))
    }else{
      return(NA)
    }
  })

  DFtofill2<-DFtofill2[!is.na(DFtofill2$Correlation),]

  if(!simple){
    DFtofill2$Comparison=factor(paste(DFtofill2$E1, "&", DFtofill2$E2))
  }else{
    sep=paste0('[',sep,']')
    DFtofill2$Comparison=factor(sapply(DFtofill2$E1,function(x) {
      tmpo<-strsplit(x,sep)[[1]][1]
      substr(tmpo,1,nchar(tmpo)-1)
    }))
  }

  if(missing(Colorpalette)){
    if(nrow(DFtofill>13)){
      autopalette<-TRUE
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
      Colorpalette<-getPalette(nrow(DFtofill))
    }
  }

  p<-ggplot(DFtofill2,aes(Cut,Correlation,group=Comparison,color=Comparison,fill=Comparison))
  if(missing(unit)){
    p<-p+xlab(XLAB)
  }else{
    p<-p+xlab(paste(XLAB,unit))
  }
  if(!newylab){
    p<-p+ylab(paste(simpleCap(method),'correlation coefficient'))
  }else{
    p<-p+ylab(ylabTitle)
  }
  p<-p+geom_rangeframe(color='black')
  if(point)  p<-p+geom_point(alpha=0.5)
  if(line)   p<-p+geom_line(alpha=0.6)
  if(smooth) p<-p+geom_smooth(method = "lm")
  if(logscaleX) p<-p + scale_x_log10()
  p<-p+theme_minimaliste(base_family = fontfamily,base_size = fontSize )
  #if(autopalette){
  p<- p + scale_color_manual(values=Colorpalette)
  p<- p + scale_fill_manual(values=Colorpalette)
  #}else{
  #  p<-p+scale_colour_brewer(palette = "Set3")
  #  p<-p+scale_fill_brewer(palette = "Set3")
  #}

  if(noTitleLegend){ p <- p + theme(legend.title = element_blank())}
  if(!missing(Title)){ p <- p + ggtitle(Title) }

  p<-p+guides(fill=guide_legend(nrow=legendN,position=position.legend),
              color=guide_legend(nrow=legendN,position=position.legend))
  p<-p+theme(legend.position = position.legend)
  ggsave(filename = filename, plot = p, device=cairo_pdf, family=fontfamily)
  print(p)
}

#' Plot the correlation of the gene expression
#' while taking in account the coefficient of variation of the genes
#' between the same tissues of two studies
#' as a function of the number of genes included in the correlation computation which
#' varies based on the minimal threshold of expression of the genes.
#'
#' @param DFa data.frame that contains the expression data for one tissue across several studies
#' @param cvDF data.frame that contains the coefficient of variation of the different genes across the tissues
#' @param sepSeq positive integer vector. Give the indices of the genes for which the expression
#'               should be plotted once ordered.
#' @param step positive integer. The step to use between the rank to plot.
#' @param method a character string indicating which correlation coefficient is to be computed.
#' One of "pearson", "kendall" or "spearman".
#' @param use character string which allows to pick the way missing data are handled.
#'            one of the strings "everything", "all.obs", "complete.obs", "na.or.complete", or "pairwise.complete.obs".
#' @param filename character string. Path where the figure should be saved.
#' @param fontfamily  character string, name of the font to use. Default: "Linux Libertine"
#' @param Colorpalette Named vector of character strings that allows to customise the colours in the plot
#' @param legend boolean. Default: FALSE. Whether the legend should be added to the plot.
#' @param point boolean. Default: FALSE. Whether to plot with geom_point()
#' @param line boolean. Default: FALSE. Whether to plot with geom_line()
#' @param smooth boolean. Default: TRUE. Whether to use a smooth line.
#'
#' @return a plot
#' @export
#'
evolCvCorrCumul<-function(DFa,cvDF,sepSeq,step=0.01,method='pearson',use='pairwise.complete.obs',
                          filename,fontfamily='Linux Libertine',Colorpalette,legend=FALSE,
                          point=FALSE,line=FALSE,smooth=TRUE){

  autopalette<-FALSE

  if(missing(sepSeq)) sepSeq<-c(seq(min(cvDF),max(cvDF),step))


  DFtofill<-as.data.frame(t(data.frame(combn(colnames(DFa),2,simplify=FALSE))),
                          stringsAsFactors = FALSE)
  colnames(DFtofill)<-c('E1','E2')
  rownames(DFtofill)<-1:nrow(DFtofill)

  DFtofill2<-DFtofill[unlist(lapply(1:nrow(DFtofill), function(x){
    if(strsplit(DFtofill[x,'E1'],split='[ (]')[[1]][1]  ==  strsplit(DFtofill[x,'E2'],split='[ (]')[[1]][1]){
      return (x)
    }else{
      return(NA)
    }
  })),]

  DFtofill2<-data.frame(na.omit(DFtofill2))
  DFtofill2<-DFtofill2[order(DFtofill2$E1),]

  if(length(sepSeq)<2){
    print("sequence of breaks too small a random value has been picked")
    sepSeq[2]<-sample(cvDF[,sample(ncol(cvDF),1)],1)
  }


  DFtofill2$Cut<-sepSeq[1]
  for(i in sepSeq[2:length(sepSeq)]){
    DFtoFillTemp<-DFtofill
    DFtoFillTemp$Cut<-i
    DFtofill2<-rbind(DFtofill2,
                     DFtoFillTemp)
  }

  e1<-e2<-NA
  DFtofill2$Correlation=sapply(1:nrow(DFtofill2),function(x){
    List[e1,e2]<-strsplit(unlist(DFtofill2[x,c('E1','E2')]), split='[ (]')
    List[e1,e2]<-lapply(list(e1,e2),function(x) return(gsub(')','',x[3])))

    #GenesID<-intersect(rownames(cvDF[cvDF[,e1] %>=% DFtofill2[x,"Cut"],]),
    #                   rownames(cvDF[cvDF[,e2] %>=% DFtofill2[x,"Cut"],]))

    averageCV<-rowMeans(cvDF[,c(e1,e2)])
    GenesID<-names(averageCV[averageCV %>=% DFtofill2[x,"Cut"]])

    temp<-DFa[GenesID,unlist(DFtofill2[x,c('E1','E2')])]

    if(nrow(temp)>0){
      return(cor(temp[,1],temp[,2],method=method,use=use))
    }else{
      return(NA)
    }
  })

  DFtofill2<-DFtofill2[!is.na(DFtofill2$Correlation),]
  print(nrow(DFtofill))
  if(missing(Colorpalette)){
    if(nrow(DFtofill>13)){
      autopalette<-TRUE
      getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Paired"))
      Colorpalette<-getPalette(nrow(DFtofill))
    }
  }

  DFtofill2$Comparison=factor(paste(DFtofill2$E1, "&", DFtofill2$E2))
  p<-ggplot(DFtofill2,aes(Cut,Correlation,group=Comparison,color=Comparison,fill=Comparison))
  p<-p+xlab('Coefficient of variation cut-off')
  p<-p+geom_rangeframe(color='black')
  if(line)  p <- p + geom_line(alpha=0.7,linetype = 4 )
  if(point) p <- p + geom_point(alpha=0.4)
  if(smooth)p <- p + geom_smooth()

  p<-p+theme_minimaliste(base_family = fontfamily)

  if (autopalette){
    p<-p+scale_color_manual(values=Colorpalette)
    p<-p+scale_fill_manual(values=Colorpalette)
  }else{
    p<-p+scale_colour_brewer(palette = "Paired")
    p<-p+scale_fill_brewer(palette = "Paired")
  }

  if(legend){
    p<-p+guides(fill=guide_legend(nrow=3,byrow=TRUE,position='bottom'))
  }else{
    p<-p+theme(legend.position="none")
  }

  print(p)
  ggsave(filename = filename, plot = p, device=cairo_pdf, family=fontfamily)
  return(p)
}




#' Plot the course of TS proteins when the level of genes correlation varies
#'
#' @param cor1 named numeric vector that comprises the correlation between two studies.
#' @param breadthExp named postive integer that comprises the breadth of expression of the proteins
#' @param simul1  named numeric vector that comprises the result of the random permutations on cor1
#' @param reference named numeric vector that comprises the correlation
#'                  between two (transcriptomic) studies that can be used as a reference
#' @param cor2 named numeric vector that comprises the correlation between two studies.
#' @param simul2 named numeric vector that comprises the result of the random permutations on cor2
#' @param aggregMethod character string for the name of the method to aggregate the different simulation together.
#'                     Default: 'rowMeans'
#' @param Title character string. Title of the plot.
#' @param CorrelationPlot boolean. Default:TRUE. Whether to add the correlation plot to help for the interpretation.
#' @param relabel  boolean. Default: TRUE. Whether the categories should be relaballed with the content of labelVec
#' @param labelVec named character string that allows to name the different categories to be displayed in the legend.
#'                 The names should be 'cor1, 'ref', 'simul1', .... and the content the labels
#'                 Default: "setNames(c("Protein/mRNA pairs",  "mRNA/mRNA pairs ('ideal' reference)",
#'                           "Randomised Protein/mRNA pairs"), c('Cor','reference','simul1'))"
#' @param sizeLine positive numeric. Size of the lines of the plot
#' @param palette named vector of character strings. Allows to customise the colours on the plot.
#'                Names should be identical to the categories or to the content of labelVec.
#' @param legendText positive numeric. Size of the legend text. Default: 0.92
#' @param base_family character string. Name of the font to use. Default: "Linux Libertine"
#' @param base_size positive numeric. Size of the font to use as a base.
#' @param centerTitle boolean. Default:TRUE. Whether to centre the title or keep it on the left.
#' @param out character string. Allows to chose the type of returned object.
#'            "plot" (default) for the plot,
#'            "shortDF" for the data.frame,
#'            "DF" for the data.frame with all the initial data for the plot,
#'            and "annotatedPlot" for an annotated plot
#'
#' @return a plot or a data.frame
#' @export
#'
cumulSpeSimulAggregated<-function(cor1,breadthExp,simul1,reference,cor2,simul2,
                                  aggregMethod='rowMeans',
                                  Title,
                                  CorrelationPlot=TRUE,
                                  relabel=TRUE,
                                  labelVec=setNames(
                                    c("Protein/mRNA pairs",
                                      "mRNA/mRNA pairs ('ideal' reference)",
                                      "Randomised Protein/mRNA pairs"),
                                    c('Cor','reference','Simul')),
                                  sizeLine=1.1,
                                  palette,legendText=0.92,
                                  base_family="Linux Libertine",
                                  base_size=11,centerTitle=TRUE,
                                  out='plot'
){

  options(stringsAsFactors=FALSE)

  #compute the fraction of specific genes for cor and for reference
  #vector already been reduced to common ids
  listDataNames<-list('cor1')
  listDataID<-list(names(cor1))
  if(!missing(reference)) {
    listDataNames[[length(listDataNames)+1]]<-'reference'
    listDataID[[length(listDataID)+1]]<-names(reference)
  }
  if(!missing(simul1)){
    listDataNames[[length(listDataNames)+1]]<-'simul1'
    listDataID[[length(listDataID)+1]]<-row.names(simul1)
  }
  if(!missing(cor2)){
    listDataNames[[length(listDataNames)+1]]<-'cor2'
    listDataID[[length(listDataID)+1]]<-names(cor2)
  }
  if(!missing(simul2)){
    listDataNames[[length(listDataNames)+1]]<-'simul2'
    listDataID[[length(listDataID)+1]]<-row.names(cor2)
  }
  listDataID[[length(listDataID)+1]]<-names(breadthExp)
  commonID<-do.call("Intersect",listDataID)
  maxSeq<-length(commonID)

  vecList<-list(simpleVecSpe(cleanupVec(cor1[commonID],commonID),breadthExp)/1:maxSeq) #collect objects to be plotted
  colNamesDF<-'Cor'

  if(!missing(simul1)){
    if(is.data.frame(simul1)){
      tempList<-lapply(colnames(simul1),function(x){
        return(simpleVecSpe(cleanupVec(setNames(simul1[commonID,x],rownames(simul1[commonID,])),commonID),breadthExp))
      })
      tempList<-as.data.frame(tempList)
      vecList[[length(vecList)+1]]<-eval(call(name=aggregMethod,tempList))/1:maxSeq
      if(!missing(simul2)){
        colNamesDF[length(colNamesDF)+1]<-"simul1"
      }else{
        colNamesDF[length(colNamesDF)+1]<-"Simul"
      }
    }else{
      vecList[[length(vecList)+1]]<-simpleVecSpe(cleanupVec(simul1,commonID),breadthExp)/1:maxSeq
      if(!missing(simul2)){
        colNamesDF[length(colNamesDF)+1]<-"simul1"
      }else{
        colNamesDF[length(colNamesDF)+1]<-"Simul"
      }
    }
  }

  if(!missing(reference)){
    vecList[[length(vecList)+1]]<-simpleVecSpe(cleanupVec(reference,commonID),breadthExp)/1:maxSeq
    colNamesDF[length(colNamesDF)+1]<-"reference"
  }

  if(!missing(cor2)){
    vecList[[length(vecList)+1]]<-simpleVecSpe(cleanupVec(cor2,commonID),breadthExp)/1:maxSeq
    colNamesDF[length(colNamesDF)+1]<-"cor2"
  }

  if(!missing(simul2)){
    if(is.data.frame(simul2)){
      tempList<-lapply(colnames(simul2),function(x){
        return(simpleVecSpe(cleanupVec(setNames(simul2[commonID,x],rownames(simul2[commonID,])),commonID),breadthExp))
      })
      vecList[[length(vecList)+1]]<-eval(call(name=aggregMethod,tempList))/1:maxSeq
      colNamesDF[length(colNamesDF)+1]<-"simul2"
    }else{
      vecList[[length(vecList)+1]]<-simpleVecSpe(cleanupVec(simul2,commonID),breadthExp)/1:maxSeq
      colNamesDF[length(colNamesDF)+1]<-"simul2"
    }
  }

  DF<-as.data.frame(vecList,stringsAsFactors=FALSE)
  colnames(DF)<-colNamesDF

  if(relabel){
    colnames(DF)<-sapply(colnames(DF), function(x) labelVec[x])
  }
  DF$pos<-rownames(DF)

  plotDF<-melt(DF,id.vars='pos')
  colnames(plotDF)<-c('Rank','Comparison','TSpercent')

  plotDF$Rank<-as.numeric(as.character(plotDF$Rank))


  p <- ggplot(plotDF,aes(x=Rank,y=TSpercent*100,color=Comparison,group=Comparison))+geom_line()
  p <- p + theme_bw(base_family=base_family,base_size=base_size)+theme(legend.position='bottom')
  p <- p + labs(title = Title, y="% of TS proteins",
                x="Nb of considered genes\n(ranked by decreasing order of their correlation coefficient)")
  if(centerTitle) p <- p +theme(plot.title=element_text(hjust=0.5))

  if(!missing(cor2)){
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(ncol=2))
  }else{
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(title.position = "left",ncol=1))
  }

  if(!missing(palette)) p <- p + scale_colour_manual(values=palette)

  if(missing(palette)){
    if(missing(cor2)){
      ColPalette<-setNames(c('plum','grey68','forestgreen'),
                           c("Protein/mRNA pairs",
                             "Randomised Protein/mRNA pairs",
                             "mRNA/mRNA pairs ('ideal' reference)"))
    }else{
      ColPalette<-setNames(c('plum','grey68','plum1','grey81','lightgreen'),
                           c(paste0('Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             "Reference: Uhl\u00e9n mRNA/GTEx mRNA pairs ('ideal' case)"))
    }

    if(all(colnames(DF)[-ncol(DF)] %in% names(ColPalette))) {
      p <- p + scale_colour_manual(values=ColPalette)
      palette<-ColPalette
    }
  }

  if(CorrelationPlot){
    cor1<-cor1[commonID]
    additionPlot<-data.frame(cor=sort(cor1[!is.na(cor1)],decreasing = TRUE),
                             rank=1:length(cor1[!is.na(cor1)]))
    p1 <- ggplot(additionPlot,aes(x=rank,y=cor))+geom_line(size=1,color='plum')+ylab(label = 'Correlation')+xlab('Ranks')
    p1 <- p1 + theme_classic(base_family=base_family,base_size=base_size)+theme(legend.position='none',plot.margin=unit(c(0.3,0,-1,0), "cm"))
    p1 <- p1 + geom_rangeframe(color='black')
    Inter05=sum(cor1>0.5)
    Inter0=sum(cor1>0)
    p1<-p1+geom_vline(colour='steelblue1',xintercept=Inter05)
    p1<-p1+geom_vline(colour='slategrey',xintercept=Inter0)
    #p1<-p1+geom_hline(colour='steelblue1',yintercept=0.5)
    #p1<-p1+geom_hline(colour='slategrey',yintercept=0)
    p1<-p1+geom_segment(y=0,yend=0,x=0,xend=length(cor1),color='slategrey',size=0.5,alpha=0.5)
    p1<-p1+geom_segment(y=0.5,yend=0.5,x=0,xend=length(cor1),color='steelblue1',size=0.5,alpha=0.5)

    if(out=='plot') print(p1)
    p<-p+geom_vline(colour='steelblue1',xintercept=Inter05)
    p<-p+geom_vline(colour='slategrey',xintercept=Inter0)
    p<-p+theme(plot.margin=unit(c(0,0,0,0), "cm"))
    #print(p)
  }
  switch(out,
         "shortDF"=return(DF),
         "DF" = return(plotDF),
         "plot"= return(p),
         "annotatedPlot"={
           gp<-suppressWarnings(ggplot_gtable(ggplot_build(p)))
           gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
           gp1$widths<-gp$widths
           return(suppressWarnings(grid.arrange(arrangeGrob(gp1, gp, heights=c(1.3,5)))))
         }
  )
}


#' Plot a figure based on the data.frame created by the cumulSpe-like function
#'
#' @param DF_long data.frame in a long format
#' @param cor1 named numeric vector that comprises the correlation between two studies.
#' @param reference named numeric vector that comprises the correlation
#'                  between two (transcriptomic) studies that can be used as a reference
#' @param Title character string. Title of the plot.
#' @param CorrelationPlot boolean. Default:TRUE. Whether to add the correlation plot to help for the interpretation.
#' @param relabel boolean. Default: TRUE. Whether the categories should be relaballed with the content of labelVec
#' @param cor2 named numeric vector that comprises the correlation between two studies.
#' @param simul1 named numeric vector that comprises the result of the random permutations on cor1
#' @param aggregateBool boolean. Default: TRUE. Whether to aggregate all the simulation of a given comparison together.
#' @param aggregMethod character string for the name of the method to aggregate the different simulation together.
#'                     Default: 'rowMeans'
#' @param labelVec named character string that allows to name the different categories to be displayed in the legend.
#'                 The names should be 'cor1, 'ref', 'simul1', .... and the content the labels
#'                 Default: "setNames(c("Protein/mRNA pairs",  "mRNA/mRNA pairs ('ideal' reference)",
#'                           "Randomised Protein/mRNA pairs"), c('Cor','reference','simul1'))"
#' @param sizeLine positive numeric. Size of the lines of the plot
#' @param palette named vector of character strings. Allows to customise the colours on the plot.
#'                Names should be identical to the categories or to the content of labelVec.
#' @param legendText positive numeric. Size of the legend text. Default: 0.92
#' @param base_family  character string. Name of the font to use. Default: "Linux Libertine"
#' @param base_size positive numeric. Size of the font to use as a base.
#' @param centerTitle boolean. Default:TRUE. Whether to centre the title or keep it on the left.
#' @param out character string. Allows to chose the type of returned object.
#'            "plot" (default) for the plot
#'            and "annotatedPlot" for an annotated plot
#' @param VinterO8 boolean. Default: FALSE. Add a vertical line for 80 percent of the data
#'
#' @return a figure
#' @export
#'
replotCumulSpeSimulAggregated<-function(DF_long,cor1,reference,Title,
                                        CorrelationPlot=TRUE,
                                        relabel=TRUE,cor2,
                                        simul1,aggregateBool=TRUE,
                                        aggregMethod='rowMeans',
                                        labelVec=setNames(
                                          c("Protein/mRNA pairs",
                                            "mRNA/mRNA pairs ('ideal' reference)",
                                            "Randomised Protein/mRNA pairs"),
                                          c('Cor','reference','Simul')),
                                        sizeLine=1,
                                        palette,legendText=0.92,
                                        base_family="Linux Libertine",
                                        base_size=11,centerTitle=TRUE,
                                        out='plot',VinterO8=FALSE){


  p <- ggplot(DF_long,aes(x=Rank,y=TSpercent*100,color=Comparison,group=Comparison))+geom_line()
  p <- p + theme_bw(base_family=base_family,base_size=base_size)+theme(legend.position='bottom')
  p <- p + labs(title = Title, y="% of TS proteins",
                x="Nb of considered genes\n(ranked by decreasing order of their correlation coefficient)")
  if(centerTitle) p <- p +theme(plot.title=element_text(hjust=0.5))

  if(!missing(cor2)){
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(ncol=2))
  }else{
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(-6,0,0,1)))+guides(color=guide_legend(title.position = "left",ncol=1))
  }

  if(!missing(palette)) p <- p + scale_colour_manual(values=palette)
  if(missing(palette)){
    if(missing(cor2)){
      ColPalette<-setNames(c('plum','grey68','forestgreen'),
                           c("Protein/mRNA pairs",
                             "Randomised Protein/mRNA pairs",
                             "mRNA/mRNA pairs ('ideal' reference)"))
    }else{
      ColPalette<-setNames(c('plum','grey68','plum1','grey81','lightgreen'),
                           c(paste0('Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             "Reference: Uhl\u00e9n mRNA/GTEx mRNA pairs ('ideal' case)"))
    }

    if(all(unique(DF_long$Comparison) %in% names(ColPalette))) {
      p <- p + scale_colour_manual(values=ColPalette)
      palette<-ColPalette
    }
  }

  if(CorrelationPlot){
    commonID=Intersect(names(cor1),names(reference))
    cor1<-cor1[commonID]
    cor1<-cor1[!is.na(cor1)]
    commonID<-names(cor1)
    index<-1:length(cor1)

    categoriesNames<-names(ColPalette)

    if(!missing(simul1)){
      CorDF<-data.frame(Ranked.cor=index,
                        Cor=sort(cor1[commonID],decreasing=TRUE),
                        Comparison=categoriesNames[1])

      if(!aggregateBool) {
        CorDF<-rbind(CorDF,
                     data.frame(Ranked.cor=index,
                                Cor=sort(simul1[commonID,sample(ncol(simul1),1)],decreasing=TRUE),
                                Comparison=categoriesNames[2],stringsAsFactors =FALSE))
      }else{
        simul1<-simul1[commonID,]
        sim1<-apply(simul1,2,sort)
        sim1<-eval(call(name=aggregMethod,as.matrix(sim1)))
        sim1<-sort(sim1,decreasing = TRUE)
        CorDF<-rbind(CorDF,
                     data.frame(Ranked.cor=index,
                                Cor=sim1,
                                Comparison=categoriesNames[2],stringsAsFactors =FALSE))
      }
      CorDF<-rbind(CorDF,
                   data.frame(Ranked.cor=index,
                              Cor=sort(reference[commonID],decreasing=TRUE),
                              Comparison=categoriesNames[length(categoriesNames)]))

      p1<-ggplot(CorDF,aes(x=as.numeric(Ranked.cor),y=Cor,group=Comparison,colour=Comparison))+geom_line(size=1)+ylab(label = 'Correlation')+xlab('Ranks')
    }else{
      additionPlot<-data.frame(cor=sort(cor1,decreasing = TRUE),
                               rank=1:length(cor1))
      p1 <- ggplot(additionPlot,aes(x=rank,y=cor))+geom_line(size=1,color='plum')+ylab(label = 'Correlation')+xlab('Ranks')
    }

    p1 <- p1 + theme_classic(base_family=base_family,base_size=base_size)+theme(legend.position='none',plot.margin=unit(c(0.3,0,-1,0), "cm"))
    p1 <- p1 + geom_rangeframe(color='black')
    p1 <- p1 + scale_colour_manual(values=palette)
    Inter05=sum(cor1>0.5)
    Inter08=sum(cor1>0.8)
    Inter0=sum(cor1>0)
    p1<-p1+geom_vline(colour='steelblue1',xintercept=Inter05)
    p1<-p1+geom_vline(colour='slategrey',xintercept=Inter0)
    if(VinterO8) p1<-p1+geom_vline(colour='royalblue4',xintercept = Inter08)
    #p1<-p1+geom_hline(colour='steelblue1',yintercept=0.5)
    #p1<-p1+geom_hline(colour='slategrey',yintercept=0)
    p1<-p1+geom_segment(y=0,yend=0,x=0,xend=length(cor1),color='slategrey',size=0.5,alpha=0.5)
    p1<-p1+geom_segment(y=0.5,yend=0.5,x=0,xend=length(cor1),color='steelblue1',size=0.5,alpha=0.5)
    p1<-p1+geom_segment(y=0.8,yend=0.8,x=0,xend=length(cor1),color='royalblue4',size=0.2,alpha=0.5)

    if(out=='plot') print(p1)
    p<-p+geom_vline(colour='steelblue1',xintercept=Inter05)
    p<-p+geom_vline(colour='slategrey',xintercept=Inter0)
    if(VinterO8) p<-p+geom_vline(colour='royalblue4',xintercept=Inter08)
    p<-p+theme(plot.margin=unit(c(-0.3,0,0,0), "cm"))
    #print(p)
  }

  switch(out,
         "plot"= return(p),
         "annotatedPlot"={
           gp<-suppressWarnings(ggplot_gtable(ggplot_build(p)))
           gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
           gp1$widths<-gp$widths
           return(suppressWarnings(grid.arrange(arrangeGrob(gp1, gp, heights=c(1.3,5)))))
         }
  )
}

#' Allows to plot the amount of TS genes among the genes that present
#' the highest observed correlations for gene expression between two studies
#' once their are ranked by decreasing order.
#' An additional plot presenting their correlations can be added to ease the interpretation
#'
#' @param cor1 named numeric vector that comprises the correlation between two studies.
#' @param protBreadth named postive integer that comprises the breadth of expression of the proteins
#' @param simul1 named numeric vector that comprises the result of the random permutations on cor1
#' @param reference named numeric vector that comprises the correlation
#'                  between two (transcriptomic) studies that can be used as a reference
#' @param cor2 named numeric vector that comprises the correlation between two studies.
#' @param simul2 named numeric vector that comprises the result of the random permutations on cor2
#' @param labelVec named character string that allows to name the different categories to be displayed in the legend.
#'                 The names should be 'cor1, 'ref', 'simul1', .... and the content the labels
#'                 Default: "setNames(c("Protein/mRNA pairs",  "mRNA/mRNA pairs ('ideal' reference)",
#'                           "Randomised Protein/mRNA pairs"), c('Cor','reference','simul1'))"
#' @param relabel boolean. Default: TRUE. Whether the categories should be relaballed with the content of labelVec
#' @param sizeLine positive numeric. Size of the lines of the plot
#' @param palette named vector of character strings. Allows to customise the colours on the plot.
#'                Names should be identical to the categories or to the content of labelVec.
#' @param legendText positive numeric. Size of the legend text. Default: 0.92
#' @param Title character string. Title of the plot.
#' @param bin positive integer.  For the additional plot that uses geom_point(), the interval of ranks to consider between two plotted genes.
#' @param seed numeric. Random seed. If simul1 and simul2 are not aggregation of all the random shuffling and one want to plot a specific occurence,
#'             it allows to plot always the same for a given seed
#' @param CorrelationPlot boolean. Default:TRUE. Whether to add the correlation plot to help for the interpretation.
#' @param intermediaryPlot boolean. Default:FALSE. Additional plot with the percent of TS protein for a given number of (ranked) genes
#' @param base_family character string. Name of the font to use. Default: "Linux Libertine"
#' @param base_size positive numeric. Size of the font to use as a base.
#' @param centerTitle boolean. Default:TRUE. Whether to centre the title or keep it on the left.
#' @param out character string. Allows to chose the type of returned object.
#'            "plot" (default) for the plot,
#'            "DF" for the data.frame with all the initial data,
#'            and "annotatedPlot" for an annotated plot
#' @return a plot or a data.frame
#' @export
#'
cumulSpe<-function(cor1,protBreadth,
                   simul1,reference,
                   cor2,simul2,
                   labelVec=setNames(
                     c("Protein/mRNA pairs",
                       "mRNA/mRNA pairs ('ideal' reference)",
                       "Randomised Protein/mRNA pairs"),
                     c('Cor','reference','simul1')),
                   relabel=TRUE, sizeLine=1.1,
                   palette,legendText=0.92,
                   Title,bin=100,seed=1323,
                   CorrelationPlot=TRUE,
                   intermediaryPlot=FALSE,
                   base_family="Linux Libertine",
                   base_size=11,centerTitle=TRUE,
                   out='plot'
){

  options(stringsAsFactors=FALSE)
  set.seed(seed)

  commonID<-names(cor1)

  if(!missing(reference)){
    commonID<-base::intersect(commonID,names(reference))

  }
  if(!missing(simul1)){
    if(is.data.frame(simul1)){
      commonID<-base::intersect(commonID,rownames(simul1))
      simul1<-setNames(simul1[commonID,sample(ncol(simul1),1)],rownames(simul1[commonID,]))
    }else{
      commonID<-base::intersect(commonID,names(simul1))
    }
  }

  if(!missing(cor2)){
    commonID<-base::intersect(commonID,names(cor2))
  }
  if(!missing(simul2)){
    if(is.data.frame(simul2)){
      commonID<-base::intersect(commonID,rownames(simul2))
      simul2<-setNames(simul2[commonID,sample(ncol(simul2),1)],rownames(simul2[commonID,]))
    }else{
      commonID<-base::intersect(commonID,names(simul2))
    }
  }

  maxSeq<-length(commonID)

  cor1<-cor1[commonID]
  vecList<-list(cor1[!is.na(cor1)])
  colNamesDF<-'Cor'
  if(!missing(simul1)){
    simul1<-simul1[commonID]
    vecList[[length(vecList)+1]]<-simul1[!is.na(simul1)]
    colNamesDF[length(colNamesDF)+1]<-"simul1"
  }
  if(!missing(reference)){
    reference<-reference[commonID]
    vecList[[length(vecList)+1]]<-reference[!is.na(reference)]
    colNamesDF[length(colNamesDF)+1]<-"reference"
  }

  if(!missing(cor2)){
    cor2<-cor2[commonID]
    vecList[[length(vecList)+1]]<-cor2[!is.na(cor2)]
    colNamesDF[length(colNamesDF)+1]<-"cor2"
  }
  if(!missing(simul2)){
    simul2<-simul2[commonID]
    vecList[[length(vecList)+1]]<-simul2[!is.na(simul2)]
    colNamesDF[length(colNamesDF)+1]<-"simul2"
  }

  DF<-as.data.frame(vecList)
  colnames(DF)<-colNamesDF

  lisDF<-lapply(colNamesDF,function(x){
    newDF<-DF[order(-DF[,x]),]
    indexGenes<-rownames(newDF)
    Spec<-protBreadth[indexGenes]
    Spec<-Spec[!is.na(Spec)]
    Bins<-1:(maxSeq %/% bin +1)*bin
    Cumul<-sapply(Bins,function(x) {
      if(x<length(Spec)) {
        bmax<-x
      }else{
        bmax<-length(Spec)
      }
      binf<-x+1-bin
      sum(Spec[binf:bmax]==1)
    })
    Cumul<-cumsum(Cumul)

    Cumulp<-(Cumul/Bins)*100
    if(intermediaryPlot){
      myData<-setNames(Cumulp,Bins)
      toSaveAsPDF<-graphics::barplot(myData,names.arg="",
                                     xlab="Genes seq",ylab="%",
                                     main=paste0('% of Specific proteins for ',x,'(',Title,')'))
      vps <- gridBase::baseViewports()
      pushViewport(vps$inner, vps$figure, vps$plot)
      grid.text(names(myData),
                x = unit(toSaveAsPDF, "native"), y=unit(-1, "lines"),
                just="right", rot=50)

      popViewport(3)
    }
    data.frame(Comparison=x,Bins=Bins,Cumul=Cumul,Cumulp=Cumulp)
  })

  plotDF<-Reduce(rbind,lisDF)


  if(relabel){
    #lapply(1:length(labelVec),function(x){
    #plotDF[plotDF$Comparison==names(labelVec[x]),]$Comparison<-labelVec[x]
    #})
    colectingLabel<-c()
    if('Cor' %in% names(labelVec)) {
      plotDF[plotDF$Comparison=='Cor',]$Comparison<-labelVec['Cor']
      colectingLabel[length(colectingLabel)+1]<-labelVec['Cor']
    }
    if('reference' %in% names(labelVec)){
      plotDF[plotDF$Comparison=='reference',]$Comparison<-labelVec['reference']
      colectingLabel[length(colectingLabel)+1]<-labelVec['reference']
    }
    if('simul1' %in% names(labelVec)){
      plotDF[plotDF$Comparison=='simul1',]$Comparison<-labelVec['simul1']
      colectingLabel[length(colectingLabel)+1]<-labelVec['simul1']
    }
    if('cor2' %in% names(labelVec)){
      plotDF[plotDF$Comparison=='cor2',]$Comparison<-labelVec['cor2']
      colectingLabel[length(colectingLabel)+1]<-labelVec['cor2']
    }
    if('simul2' %in% names(labelVec)){
      plotDF[plotDF$Comparison=='simul2',]$Comparison<-labelVec['simul2']
      colectingLabel[length(colectingLabel)+1]<-labelVec['simul2']
    }
  }

  p <- ggplot(as.data.frame(plotDF),aes(x=Bins,y=Cumulp,color=Comparison))+geom_line(size=sizeLine)
  p <- p + theme_bw(base_family=base_family,base_size=base_size)+theme(legend.position='bottom')
  p <- p + labs(title = Title, y="% of TS proteins",
                x="Nb of considered genes\n(ranked by decreasing order of their correlation coefficient)")
  if(centerTitle) p <- p +theme(plot.title=element_text(hjust=0.5))

  if(!missing(cor2)){
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(ncol=2))
  }else{
    p <- p + theme(legend.text= element_text(size = rel(legendText)),
                   legend.margin=margin(c(0,0,0,1)))+guides(color=guide_legend(title.position = "left",ncol=1))
  }

  if(!missing(palette)) p <- p + scale_colour_manual(values=palette)

  if(missing(palette)){
    if(missing(cor2)){
      ColPalette<-setNames(c('plum','grey68','forestgreen'),
                           c("Protein/mRNA pairs",
                             "Randomised Protein/mRNA pairs",
                             "mRNA/mRNA pairs ('ideal' reference)"))
    }else{
      ColPalette<-setNames(c('plum','grey68','plum1','grey81','lightgreen'),
                           c(paste0('Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/Uhl\u00e9n'),
                             paste0('Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             paste0('Randomised Protein/mRNA pairs: Pandey',unit,'/GTEx'),
                             "Reference: Uhl\u00e9n mRNA/GTEx mRNA pairs ('ideal' case)"))
    }
    if(all(colectingLabel %in% names(ColPalette))) {
      p <- p + scale_colour_manual(values=ColPalette)
      palette<-ColPalette
    }
  }

  if(CorrelationPlot){
    cor1<-cor1[commonID]
    additionPlot<-data.frame(cor=sort(cor1[!is.na(cor1)],decreasing = TRUE),
                             rank=1:length(cor1[!is.na(cor1)]))
    p1 <- ggplot(additionPlot,aes(x=rank,y=cor))+geom_line(size=1,color='plum')+ylab(label = 'Correlation')+xlab('Ranks')
    p1 <- p1 + theme_classic(base_family=base_family,base_size=base_size)+theme(legend.position='none',plot.margin=unit(c(0.3,0,-1,0), "cm"))
    p1 <- p1 + geom_rangeframe(color='black')

    Inter05=sum(cor1>0.5)
    Inter0=sum(cor1>0)
    p1<-p1+geom_vline(colour='steelblue1',xintercept=Inter05)
    p1<-p1+geom_vline(colour='slategrey',xintercept=Inter0)
    #p1<-p1+geom_hline(colour='steelblue1',yintercept=0.5)
    #p1<-p1+geom_hline(colour='slategrey',yintercept=0)
    p1<-p1+geom_segment(y=0,yend=0,x=0,xend=length(cor1),color='slategrey',size=0.5,alpha=0.5)
    p1<-p1+geom_segment(y=0.5,yend=0.5,x=0,xend=length(cor1),color='steelblue1',size=0.5,alpha=0.5)

    if(out=='plot') print(p1)
    p<-p+geom_vline(colour='steelblue1',xintercept=Inter05)
    p<-p+geom_vline(colour='slategrey',xintercept=Inter0)
    p<-p+theme(plot.margin=unit(c(0,0,0,0), "cm"))
  }

  #print(p)
  switch (out,
          "DF" = return(plotDF),
          "plot"= return(p),
          "annotatedPlot"={
            gp<-suppressWarnings(ggplot_gtable(ggplot_build(p)))
            gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
            gp1$widths<-gp$widths
            return(suppressWarnings(grid.arrange(arrangeGrob(gp1, gp, heights=c(1.3,5)))))
          }
  )
}


#' Rank the genes based on their decreasing order of correlation and then plot them
#'
#' @param DF.raw numeric data.frame that contains all the correlation
#' @param Title character string. Title to figure on the plot.
#' @param bin positive integer. eg if N, 1 point out of N points is plotted.
#' @param Table data.frame that compiles many information
#'
#' @return a plot
#' @export
#'
Cumulative_spe<-function(DF.raw,Title,bin=100,Table){

  options(stringsAsFactors=FALSE)

  maxSeq<-nrow(DF.raw)
  DF<-DF.raw[,-ncol(DF.raw)]
  #index<-1:maxSeq
  colNamesDF<-colnames(DF)

  lisDF<-lapply(colNamesDF,function(x){
    newDF<-DF[order(-DF[,x]),]
    indexGenes<-rownames(newDF)
    Spec<-setNames(Table[indexGenes,]$nb.Tissues.proteomics,indexGenes)
    Bins<-1:(maxSeq%/%bin +1)*bin
    Cumul<-sapply(Bins,function(x) {
      if(x< length(Spec)) {
        bmax<-x
      }else{
        bmax<-length(Spec)
      }
      binf<-x+1-bin

      sum(Spec[binf:bmax]==1)
    })
    Cumul<-cumsum(Cumul)
    toto<-Bins
    toto[length(toto)]<-nrow(newDF)
    Cumulp<-(Cumul/toto)*100

    myData<-setNames(Cumulp,Bins)
    cairo_pdf(paste0(Title,'_',x,'.pdf'))
    toSaveAsPDF<-graphics::barplot(myData,names.arg="",
                                   xlab="Genes seq",ylab="%",
                                   main=paste0('% of Specific proteins for ',x,'(',Title,')'))
    vps <- gridBase::baseViewports()
    pushViewport(vps$inner, vps$figure, vps$plot)


    grid.text(names(myData),
              x = unit(toSaveAsPDF, "native"), y=unit(-1, "lines"),
              just="right", rot=50)

    popViewport(3)
    dev.off()


    data.frame(Comparison=x,Bins=Bins,Cumul=Cumul,Cumulp=Cumulp)
  })
  plotDF<-Reduce(rbind,lisDF)

  p<-ggplot(as.data.frame(plotDF),aes(x=Bins,y=Cumulp,color=Comparison))+geom_line()
  p<-p+theme_bw()+theme(legend.position='bottom')+ labs(title = Title)
  print(p)

  return(plotDF)
}



# Focused on the specific genes ------------------------------


#' Draw the heatmap based on the output of matrix overlap_res
#'
#' @param dfPercent overlap matrix computed by matrix overlap_res
#' @param dfpvalue associated p-value matrix
#' @param label.df1 character string. Label for the first data.frame
#' @param label.df2 character string. Label for the second data.frame
#' @param col colour palette for the main heatmap
#' @param colp colour palette for the p-value associated heatmap.
#'
#' @return two heatmaps one with the percentage of overlap between the two studies
#'         and the second with the corresponding p-values.
#' @export
#'
drawOverlapHeatmaps<-function(dfPercent,dfpvalue,label.df1,label.df2,
                              col=barzinePhdR::colCbc(),colp){

  if(missing(colp)) colp=col

  cat('\n###### Ratio\n')
  graphics::par(oma=c(2.5,2.5,2.5,2.5))
  heatmap.2(as.matrix(dfPercent), trace="none", dendrogram="none",Rowv=FALSE,Colv=FALSE,
            cellnote=signif(as.matrix(dfPercent),digits=2),
            cexRow =1,cexCol = 1, notecol="black", col=col,xlab=label.df2,ylab=label.df1)

  cat('\n\n###### p-value\n')

  graphics::par(oma=c(2.5,2.5,2.5,2.5))
  heatmap.2(as.matrix(dfpvalue), trace="none", dendrogram="none",Rowv=FALSE,Colv=FALSE,
            cellnote=signif(as.matrix(dfpvalue),digits=2),
            cexRow =1,cexCol = 1, notecol="black", col=colp, xlab=label.df2,ylab=label.df1)
}



#' Calculate the overlap of specific genes (uniquely expressed in one tissue)
#' between two studies.
#'
#' @param DF1 Expression data.frame of the first studie to be compared
#' @param DF2 Expression data.frame of the second studie to be compared
#' @param ratio numeric; multiplicator for the number of genes used
#'              for the second data.frame compared to the first one.
#'              "-1" if all specific genes of each data.frame should be compared
#' @param common argument passed on to overlapProtTrans.Specific; default: FALSE
#' @param verbose boolean; default: TRUE. Allows outputting additional information
#' @param selectSpe boolean; default: FALSE. Allows to select the specific genes in the data.frames
#'                 if not done beforehand.
#' @param report character string. Default: 'html'. Other option: ""
#'               When report='html' use cat.html instead of print
#'
#' @return a list of two matrices. One matrix with the overlap (in percentage)
#'         and another matrix with the corresponding p-values.
#' @export
#'
matrix.overlap_res<-function(DF1,DF2,ratio=-1,common=FALSE,
                             selectSpe=FALSE,verbose=TRUE,
                             report='html'){

  if(selectSpe){
    if(verbose) print('Selection of the specific genes for each data.frame\nNew nb of rows for DF1 and DF2')
    DF1<-selectSpecific(DF1,verbose=verbose)
    DF2<-selectSpecific(DF2,verbose=verbose)
  }

  if(verbose) cat('\n###### Ratio\n')
  df1_df2_Percent<-as.data.frame(Reduce(cbind,
                                        lapply(colnames(DF2),
                                               function(x){
                                                 sapply(colnames(DF1),function (y){
                                                   overlapProtTrans.Specific(DF1,DF2,y,x,ratio=ratio,
                                                                             fig=FALSE,out='percent',
                                                                             common=common,report=report,
                                                                             verbose=verbose)
                                                 })})))
  colnames(df1_df2_Percent)<-colnames(DF2)

  if(verbose) cat('\n\n###### p-value\n')
  df1_df2_pvalue<-as.data.frame(Reduce(cbind,
                                       lapply(colnames(DF2),
                                              function(x){
                                                sapply(colnames(DF1), function (y){
                                                  overlapProtTrans.Specific(DF1,DF2,y,x,ratio=ratio,
                                                                            fig=FALSE,out='p-value',
                                                                            common=common,report=report,
                                                                            verbose=verbose)
                                                })})))
  colnames(df1_df2_pvalue)<-colnames(DF2)

  return(list(df1_df2_Percent,df1_df2_pvalue))
}




#' Compute test of significance and draw venn diagram
#' of the observed overlap between the two data.frames.
#'
#' @param DF1 numeric data.frame of the first study to compare
#' @param DF2 numeric data.frame of the second study to compare
#' @param cond1 character string, column name of the first data.frame being considered for the comparison
#' @param cond2 character string, column name of the second data.frame that being considered for the comparison
#' @param ratio integer. Default: 1. Possible multiplier to apply to the first data.frame genes number
#'              to get the second data.frame genes number to consider for the comparison
#' @param common boolean. Default: TRUE.
#'               Whether the two data.frames should comprise only the same genes for the comparison
#' @param fig boolean. Default: TRUE. Whether the figure should also be printed directly
#' @param out character string. "percent" or "pvalue"
#' @param threshold numeric. Default: 1. Minimal expression for DF2 to be considered for the comparison
#' @param thresholdDF1 numeric. Default:0. Minimal expression for DF1 to be considered for the comparison
#' @param categories character string vector of two. default: c('DF1','DF2'). Allows to give the name of the studies.
#' @param Col numerical vectors of two. Default: "c('coral4','darkorchid1')" Colours for the venn diagram
#' @param report character string. Default: ''. When report='html' use cat.html instead of print
#' @param verbose boolean. Default: TRUE
#'
#' @return significance test and venn diagram
#' @export
#'
overlapProtTrans.Specific<-function(DF1,DF2,cond1,cond2,ratio=1,common=TRUE,fig=TRUE,out,
                                    threshold=1,thresholdDF1=0,categories=c('DF1','DF2'),
                                    Col=c('coral4','darkorchid1'),report='',verbose=TRUE){

  grid::grid.newpage()
  if(missing(cond2)) cond2<-cond1

  if(common){
    comNames<-base::intersect(rownames(DF1),rownames(DF2))
    DF1<-DF1[comNames,]
    DF2<-DF2[comNames,]
  }

  vecDF2.raw<-setNames(DF2[base::order(DF2[,cond2],decreasing=TRUE),][,cond2],
                       rownames(DF2[base::order(DF2[,cond2],decreasing=TRUE),]))
  vecDF1.raw<-setNames(DF1[base::order(DF1[,cond1],decreasing=TRUE),][,cond1],
                       rownames(DF1[base::order(DF1[,cond1],decreasing=TRUE),]))
  if(thresholdDF1==0){
    vecDF1<-vecDF1.raw[vecDF1.raw > thresholdDF1]
  }else{
    vecDF1<-vecDF1.raw[vecDF1.raw >= thresholdDF1]
  }

  if(verbose){
      if(report=='html'){
        cat.html(paste(cond1,cond2,sep='~'))
        cat.html(paste('Number of features tested from first data.frame given as input ',length(vecDF1)))
        }else{
          print(paste(cond1,cond2,sep='~'))
          print(paste('Number of features tested from first data.frame given as input ',length(vecDF1)))
        }
  }

  if(ratio!=-1){
    lg<-length(vecDF1)*ratio
    vecDF2<-vecDF2.raw[1:lg]
  }else{
    vecDF2<-vecDF2.raw[vecDF2.raw>=threshold]
    lg<-length(vecDF2.raw)
  }

  if(report=='html'){
    if(verbose){
      cat.html(paste('Number of features from the second data.frame ',length(vecDF2)))
      cat.html(paste('Total number of features in the second data.frame ',length(vecDF2.raw)))
    }

    cat.html('hypergeometric test')
    try(cat.html(print(phyper(length(intersect(names(vecDF2),names(vecDF1)))-1,length(vecDF1),
                              length(vecDF2.raw)-length(vecDF1),
                              length(vecDF2),lower.tail=FALSE))))
    pvalue<-2
    try(pvalue<-phyper(length(intersect(names(vecDF2),names(vecDF1)))-1,length(vecDF1),
                       length(vecDF2.raw)-length(vecDF1),
                       length(vecDF2),lower.tail=FALSE))

    cat.html('binomial (greater than)')
    try(cat.html(print(1 - pbinom(length(intersect(names(vecDF2),names(vecDF1))),
                                  lg,length(vecDF1)/length(vecDF2.raw),lower.tail=FALSE))))

    cat.html('binomial test')
    try(cat.html(print(binom.test(length(intersect(names(vecDF2),names(vecDF1))),
                                  length(vecDF2),length(vecDF1)/length(vecDF1.raw),alternative='greater'))))
  }else{
    if(verbose){
      print(paste('Number of features from the second data.frame ',length(vecDF2)))
      print(paste('Total number of features in the second data.frame ',length(vecDF2.raw)))
    }

    print('hypergeometric test')
    try(print(phyper(length(intersect(names(vecDF2),names(vecDF1)))-1,length(vecDF1),
                     length(vecDF2.raw)-length(vecDF1),
                     length(vecDF2),lower.tail=FALSE)))
    pvalue<-2
    try(pvalue<-phyper(length(intersect(names(vecDF2),names(vecDF1)))-1,length(vecDF1),
                       length(vecDF2.raw)-length(vecDF1),
                       length(vecDF2),lower.tail=FALSE))

    print('binomial (greater than)')

    try(print(1 - pbinom(length(intersect(names(vecDF2),names(vecDF1))),
                         lg,length(vecDF1)/length(vecDF2.raw),lower.tail=FALSE)))

    try(print(1 - pbinom(length(intersect(names(vecDF2),names(vecDF1))),
                         lg,length(vecDF1)/length(vecDF2.raw),lower.tail=FALSE)))

    print('binomial test')
    try(print(binom.test(length(intersect(names(vecDF2),names(vecDF1))),
                         length(vecDF2),length(vecDF1)/length(vecDF1.raw),alternative='greater')))
  }
  if(!fig){

    if(!missing(out)){
      if(out=='percent'){
        return(length(intersect(names(vecDF2),
                                names(vecDF1)))/(length(vecDF2)+length(vecDF1)-length(intersect(names(vecDF2),names(vecDF1)))))
      }
      if(out=="p-value"){
        return(pvalue)
      }
    }
    return()
  }

  pvalue<-paste('p-value = ',pvalue)
  p<-VennDiagram::draw.pairwise.venn(ind=FALSE,area1=length(vecDF2),area2=length(vecDF1),
                                     cross.area=length(intersect(names(vecDF2),names(vecDF1))),
                                     category = c(categories[2],categories[1]), col = c(Col[2],Col[1]),
                                     fill = c(Col[2],Col[1]), cat.col= c(Col[2],Col[1]),
                                     cex=2, alpha = rep(0.4, 2),euler.d=TRUE,margin=0, main=cond1)
  gridExtra::grid.arrange(grobs=list(sub=textGrob(cond1,gp=gpar(fontsize=30)),
                                     textGrob(pvalue,gp=gpar(fontsize=10)),gTree(children=p)),
                          heights=c(1,0.5,5))

}


# Analyses based on breadth of expression -----------------

#' For easier comparison, plot the uniquely expressed genes (colored by tissues) in two studies
#'
#' @param DF1 numeric data.frame (first study expression data)
#' @param DF2 numeric data.frame (second study expression data)
#' @param threshold1 numeric. Expression above which a gene is considered as expressed for the first study
#' @param threshold2 numeric. Expression above which a gene is considered as expressed for the second study
#' @param label1 character string. Label for the first study to use on the plot
#' @param label2 character string. Label for the second study to use on the plot
#' @param sorted boolean. Default: TRUE. Whether the tissues should be sorted in function of their number of tissue specific genes
#' @param common boolean. Default: TRUE. Whether the two studies should share identical rownames and colnames
#' @param colorpal colour palette to use in the figure (done with ggplot2::scale_fill_manual)
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param output character string. Switch that allows to choose between 'count' for the count of unique genes across the tissues
#'               or a ratio based on the distribution of the tissue specific genes across each study.
#' @param verbose boolean. Default: TRUE.
#' @param ... other arguments that can be used by ggplot2::theme_bw()
#'
#' @return a figure
#' @export
#'
bibarplotsDiversityCond<-function(DF1,DF2,threshold1=0,threshold2=1,
                                  label1='Proteomics (detected)',
                                  label2=paste('Transcriptomics (\u2265 ',threshold2,'FPKM)'),
                                  sorted=TRUE, common=TRUE,
                                  colorpal=NULL,publish=TRUE,
                                  output='count',
                                  verbose=TRUE,...){

  #note \u2265 = "≥"

  col1<-colnames(DF1)
  col2<-colnames(DF2)

  if(common) {
    com<-intersect(col1,col2)
    DF1<-DF1[,com]
    DF2<-DF2[,com]
  }

  d1.spe<-prep(DF1,threshold1,label1)
  d2.spe<-prep(DF2,threshold2,label2)

  yli<-max(nrow(d1.spe),nrow(d2.spe))

  if(sorted){
    d1.uvar<-unique(as.character(d1.spe$variable))
    d2.uvar<-unique(as.character(d2.spe$variable))
    nb.cond<-length(intersect(d1.uvar,d2.uvar))
    if(verbose) print(paste(nb.cond,'is the number of conditions presenting unique features in BOTH datasets'))
    sortD1<-prep2(d1.uvar,d1.spe,decreas = TRUE)
    sortD2<-prep2(d2.uvar,d2.spe,decreas = TRUE)
  }

  DF_bars<-as.data.frame(rbind(d1.spe,d2.spe))
  colnames(DF_bars)[2]<-'Tissue'
  DF_bars$rank<-0

  switch(output,
         'count'={
           toptitle='Number of unique features per tissue'
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD1))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD1))
           }
           p1 <- ggplot(DF_bars[DF_bars$Label==label1,],
                        aes(x=Label))+geom_bar(aes(y=..count..,
                                                   fill=Tissue),width=0.5)
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD2))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD2))
           }
           p2 <- ggplot(DF_bars[DF_bars$Label==label2,],
                        aes(x=Label))+geom_bar(aes(y=..count..,
                                                   fill=Tissue),width=0.5)
           yla<-'count'
         },
         {
           toptitle="Distribution of unique features per tissue"
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD1))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD1))
           }
           p1 <- ggplot(DF_bars[DF_bars$Label==label1,],
                        aes(x=Label,fill=Tissue))+geom_bar(position='fill',width=0.5)
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD2))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD2))
           }
           p2 <- ggplot(DF_bars[DF_bars$Label==label2,],
                        aes(x=Label,fill=Tissue))+geom_bar(position='fill',width=0.5)
           yla<-'ratio'
           # yli<-1
         }
  )

  p1 <- p1 + coord_flip()
  p2 <- p2 + coord_flip()

  if(!missing(colorpal)){
    p1 <- p1 + scale_fill_manual(values=colorpal)
    p2 <- p2 + scale_fill_manual(values=colorpal)
  }
  if(publish) {
    p1 <- p1 + theme_bw(...)
    p2 <- p2 + theme_bw(...)
  }
  legend<-g_legend(p1)

  p1 <- p1 + xlab(label1)+theme(legend.position='none',axis.text.y=element_blank())
  p2 <- p2 + xlab(label2)+theme(legend.position='none',axis.text.y=element_blank())

  p1 <- p1 + ylab(yla)
  p2 <- p2 + ylab(yla)

  if(output=='count'){
    p1 <- p1 + ylim(0,yli)
    p2 <- p2 + ylim(0,yli)
  }

  gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
  gp2<-suppressWarnings(ggplot_gtable(ggplot_build(p2)))

  gp2$widths <- gp1$widths
  gp2$heights <-gp1$heights

  lay<-rbind(c(1,1,1,3),c(2,2,2,3))

  grid.arrange(p1,p2,legend,layout_matrix=lay,
               top=toptitle)
}


#' For easier comparison, plot the uniquely expressed genes (colored by tissues) in two studies
#' Allows to consider several thresholds of expression for the second study.
#'
#' @param DF1 numeric data.frame (first study expression data)
#' @param DF2 numeric data.frame (second study expression data)
#' @param threshold1 numeric. Expression above which a gene is considered as expressed for the first study
#' @param threshold2 a vector of numeric. Expression values above which a gene is considered as expressed for the second study
#' @param label1 character string. Label for the first study to use on the plot (first part)
#'               Default: 'Proteins'
#' @param midLabel1 character string. Label for the first study to use on the plot (first part before threshold)
#'                  Default: '('
#' @param endLabel1 character string. Label for the first study to use on the plot (second part after threshold)
#'                  Default: ')'
#' @param label2 character string. Label for the second study to use on the plot (first part)
#'               Default: 'mRNA'
#' @param midLabel2 character string. Label for the second study to use on the plot (first part before threshold)
#'                  Default: '('
#' @param endLabel2 character string. Label for the second study to use on the plot (second part after threshold)
#'                  Default: 'FPKM)'
#' @param sorted boolean. Default: TRUE. Whether the tissues should be sorted in function of their number of tissue specific genes
#' @param common boolean. Default: TRUE. Whether the two studies should share identical rownames and colnames
#' @param colorpal colour palette to use in the figure (done with ggplot2::scale_fill_manual)
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param output character string. Switch that allows to choose between 'count' for the count of unique genes across the tissues
#'               or a ratio based on the distribution of the tissue specific genes across each study.
#' @param verbose  boolean. Default: TRUE.
#' @param decreasing boolean. Default: TRUE. Whether the sorting should be done in decreasing order.
#' @param ... other arguments that can be used by ggplot2::theme_bw()
#'
#' @return a list of object that would allow to create a figure
#' @export
#'
multibarplotsDiversityCond<-function(DF1,DF2,threshold1=0,threshold2=c(0,1,5),
                                     label1='Proteins',
                                     midLabel1='(',
                                     endLabel1=')',
                                     label2='mRNA',
                                     midLabel2='(',
                                     endLabel2='FPKM)',
                                     sorted=TRUE, common=TRUE,
                                     colorpal=NULL,publish=TRUE,
                                     output='count',
                                     verbose=TRUE,
                                     decreasing=TRUE,...){


  col1<-colnames(DF1)
  col2<-colnames(DF2)

  if(common) {
    com<-intersect(col1,col2)
    DF1<-DF1[,com]
    DF2<-DF2[,com]
  }

  d1.spe<-d2.spe<-NA

  List[d1.spe,label1]<-prep(DF1,threshold=threshold1,label=label1,midLabel=midLabel1,endLabel=endLabel1)
  d2.spe<-lapply(threshold2,function(x) {
    prep(DF2,threshold=x,label=label2,midLabel=midLabel2,endLabel=endLabel2)
  })

  label2<-lapply(d2.spe,function(x) return(x[[2]]))
  d2.spe<-lapply(d2.spe,function(x) return(x[[1]]))

  yli<-max(nrow(d1.spe),sapply(d2.spe,nrow))

  if(sorted){
    d1.uvar<-unique(as.character(d1.spe$variable))
    sortD1<-prep2(d1.uvar,d1.spe)
    sortD2<-sapply(d2.spe,function(x){
      d2.uvar<-unique(as.character(x$variable))
      nb.cond<-length(intersect(d1.uvar,d2.uvar))
      if(verbose) print(paste(nb.cond,'is the number of conditions presenting unique features in BOTH datasets'))
      return(prep2(d2.uvar,x))
    })}

  DF_bars<-data.table::rbindlist(d2.spe)
  DF_bars<-as.data.frame(rbind(d1.spe,DF_bars))
  colnames(DF_bars)[2]<-'Tissue'
  DF_bars$rank<-0

  switch(output,
         'count'={
           toptitle='Number of unique features per tissue'
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD1))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD1))
           }
           p1 <- ggplot(DF_bars[DF_bars$Label==label1,],
                        aes(x=Label))+geom_bar(aes(y=..count..,
                                                   fill=Tissue),width=0.5)
           p2.list<-lapply(1:length(d2.spe),function(y){
             if(sorted){
               DF_bars$rank<-sapply(DF_bars$Tissue,
                                    function(x){grep(x,names(sortD2[[y]]))})
               DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD2[[y]]))
             }
             p2 <- ggplot(DF_bars[DF_bars$Label==label2[[y]],],
                          aes(x=Label))+geom_bar(aes(y=..count..,
                                                     fill=Tissue),width=0.5)
             return(p2)
           })
           yla<-'count'

         },
         {
           toptitle="Distribution of unique features per tissue"
           if(sorted){
             DF_bars$rank<-sapply(DF_bars$Tissue,
                                  function(x){grep(x,names(sortD1))})
             DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD1))
           }
           p1 <- ggplot(DF_bars[DF_bars$Label==label1,],
                        aes(x=Label,fill=Tissue))+geom_bar(position='fill',width=0.5)

           p2.list<-lapply(1:length(d2.spe),function(y){
             if(sorted){
               DF_bars$rank<-sapply(DF_bars$Tissue,
                                    function(x){grep(x,names(sortD2[[y]]))})
               DF_bars$Tissue<-factor(DF_bars$Tissue,levels=names(sortD2[[y]]))
             }
             p2 <- ggplot(DF_bars[DF_bars$Label==label2[[y]],],
                          aes(x=Label,fill=Tissue))+geom_bar(position='fill',width=0.5)
             return(p2)
           })
           yla<-'ratio'
           # yli<-1
         })

  p1 <- p1 + coord_flip()
  p2.list <-lapply(p2.list, function(p3) return(p3 + coord_flip()))

  if(!missing(colorpal)){
    p1 <- p1 + scale_fill_manual(values=colorpal)
    p2.list<- lapply(p2.list, function(p3) return(p3 + scale_fill_manual(values=colorpal)))
  }
  if(publish) {
    p1 <- p1 + theme_bw(...)
    p2.list <- lapply(p2.list, function(p3) return(p3+ theme_bw(...)))
  }
  legend<-g_legend(p1)

  p1 <- p1 + xlab(label1)+theme(legend.position='none',axis.text.y=element_blank())
  p2.list <- lapply(1:length(d2.spe), function(y) return(p2.list[[y]] + xlab(label2[[y]])+theme(legend.position='none',axis.text.y=element_blank())))

  p1 <- p1 + ylab(yla)
  p2.list <- lapply(p2.list, function(p3) return(p3 + ylab(yla)))

  if(output=='count'){
    p1 <- p1 + ylim(0,yli)
    p2.list <- lapply(p2.list, function(p3) return(p3 + ylim(0,yli)))
  }

  gp1<-suppressWarnings(ggplot_gtable(ggplot_build(p1)))
  gp2.list<-lapply(p2.list,function(x) suppressWarnings(ggplot_gtable(ggplot_build(x))))

  gp2.list<-lapply(gp2.list, function(x){
    x$widths <- gp1$widths
    return(x)})
  gp2.list<-lapply(gp2.list, function(x){
    x$heights <-gp1$heights
    return(x)})

  #  lay<-rbind(c(1,1,1,3),c(2,2,2,3))


  #  grid.arrange(p1,p2,legend,layout_matrix=lay,
  #               top=toptitle)

  return(list(legend,gp1,gp2.list))
}


#' Plot the breadth of a first study compared to a second one.
#' @description The colours refer to the breadth of expression in the other study.
#'
#' @param a numeric data.frame comprising the expression data of the first study.
#' @param b numeric data.frame comprising the expression data of the second study.
#' @param lapse positive integer.
#'              Allow to relax constraints on perfect equality of breadth for the genes between the two studies.
#'              A "similar" class is created. The lapse allows to define which level of similarity is acceptable.
#'              e.g. if lapse=3, all genes that have a breadth of expression that varies at most of 3 are considered
#'              to have similar breadth.
#' @param typeR character string. Allows to pick the type of output.
#'              "vec" for a vector of the type of breath of expression observed in the first study compared to the second one.
#' @param simplify boolean. Whether to focus only on the most extreme breadth of expression.
#'                 "Expressed in all" or "nearly in all" might be more descriptive than "expressed in 5" versus "expressed in 6".
#' @param strict boolean. Whether if only strict equality between the breadth of expression of each gene between the two studies should be considered
#'               or if the "lapse" argument should be used to define 'similar' cases as well.
#' @param Lims a length-4 numeric vector. When the analysis is in the "simplify" mode,
#'             it allows to delimit two ranges of expression breadth to compare between the two studies.
#' @param colName character string. Gives the name of the column in the data.frame that has recorded the breath of expression of the genes
#'
#' @return depending on "typeR" either a vector or a data.frame
#' @export
#'
sharedBreadth<-function(a,b,lapse,typeR='vec',simplify,strict,
                        Lims, colName='nb.tissues'){

  if(missing('Lims')) Lims=c(1,3,10,12)
  a$shared<-unlist(lapply(rownames(a),function(x) {
    if(!is.na(b[x,colName])){
      if(!simplify){
        if(a[x,colName]==b[x,colName]){
          shared='Identical'
        }else{
          if(!strict){
            if(lapse >= abs(a[x,colName]-b[x,colName])){
              shared='Similar'
            }else{
              shared='Different'
            }
          }
        }
      }else{ #should be simplified
        if(any(a[x,colName]==c(Lims[1]:Lims[2],Lims[3]:Lims[4]))){
          if(a[x,colName]==b[x,colName]){
            shared='Identical'
          }else{
            if(!strict){
              if(lapse >= abs(a[x,colName]-b[x,colName])){
                shared='Similar'
              }else{
                shared='Different'
              }
            }
          }
        }else{
          shared='Mixed'
        }
      }
    }else{
      shared='Unshared'}
  }
  ))

  levelOrder<-factor(c('Identical','Similar','Mixed','Different','Unshared'))
  a$shared<-factor(a$shared,levels=levelOrder)
  switch(typeR,
         'vec'= return(a$shared),
         'df' = return(a)
  )
}


#' Compute and compare the breath of expression of two studies and can plot them
#'
#' @param DF1 numeric data.frame comprising the expression data of the first study.
#' @param DF2 numeric data.frame comprising the expression data of the second study.
#' @param omit.zero boolean. Default: TRUE.
#'                  Whether the genes which do not reach the minimal expression threshold in any tissue
#'                 should be stripped from the output.
#' @param threshold numeric. Default: 0. Minimal level of expression to be considered as expressed.
#' @param strict boolean. Default:FALSE. Whether if only strict equality between the breadth of expression of each gene between the two studies should be considered
#'               or if the "lapse" argument should be used to define 'similar' cases as well.
#' @param simplify boolean. Default: TRUE. Whether to focus only on the most extreme breadth of expression.
#'                 "Expressed in all" or "nearly in all" might be more descriptive than "expressed in 5" versus "expressed in 6".
#' @param Lims a length-4 numeric vector. Default: c(1,3,10,12).
#'             When the analysis is in the "simplify" mode,
#'             it allows to delimit two ranges of expression breadth to compare between the two studies.
#' @param lapse positive integer. Default: 3. Allow to relax constraints on perfect equality of breadth for the genes between the two studies.
#'              A "similar" class is created. The lapse allows to define which level of similarity is acceptable.
#'              e.g. if lapse=3, all genes that have a breadth of expression that varies at most of 3 are considered
#'              to have similar breadth.
#' @param colour colour palette for the different categories.
#'               Default: "set2" of colorbrewer2
#' @param ... other arguments that can be used by ggplot2::theme_bw()
#' @param annotate boolean. Default: TRUE. Whether to annotate the figures with the counts in each category for an accurate read.
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param P1 boolean. Default: TRUE. Whether to plot the figure for the first study and to return it as a port of the result.
#' @param P2 boolean. Default: TRUE. Whether to plot the figure for the second study and to return it as a port of the result.
#' @param unit.DF1 character string. Default: '' For the labelling, allows to input the correct unit in which the genes in the first study are expressed.
#' @param unit.DF2 character string. Default: '' For the labelling, allows to input the correct unit in which the genes in the second study are expressed.
#' @param Prefix character string. Default: ''. Allows to complete the labelling on the plot.
#'               e.g. Prefix='Expression'
#'
#' @return depending of whether P1 and P2 are true,
#'         the output is the data.frames with an added column for the expression breadth
#'         and the plot of their distribution.
#' @export
#'
cross.sharedDistrib<-function(DF1,DF2,omit.zero=TRUE,threshold=0,strict=FALSE,
                              simplify=TRUE,Lims=c(1,3,10,12),lapse=3,colour,...,
                              annotate=TRUE,publish=TRUE,P1=TRUE,P2=FALSE,
                              unit.DF1='',unit.DF2='',Prefix=''){
  identicalRowNames=FALSE
  if(nrow(DF1)==nrow(DF2)){
    if(all(sort(rownames(DF1)==sort(rownames(DF2))))) identicalRowNames=TRUE
  }
  if(length(threshold)==1){
    DF1<-computeBreadth(DF1,threshold=threshold,omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold,omit.zero=omit.zero,typeR='df')
  }else{
    DF1<-computeBreadth(DF1,threshold=threshold[1],omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold[2],omit.zero=omit.zero,typeR='df')
  }

  if (missing(lapse)) lapse=3
  if(!simplify){
    DF1$Genes<-sharedBreadth(DF1,DF2,simplify=simplify,strict=strict,lapse=lapse)
    DF2$Genes<-sharedBreadth(DF2,DF1,simplify=simplify,strict=strict,lapse=lapse)
  }else{
    if(missing(Lims)) Lims=c(1,3,10,12)
    DF1$Genes<-sharedBreadth(DF1,DF2,simplify=simplify,strict=strict,Lims=Lims,lapse=lapse)
    DF2$Genes<-sharedBreadth(DF2,DF1,simplify=simplify,strict=strict,Lims=Lims,lapse=lapse)
  }

  if(identicalRowNames) {
    DF1$Genes<-as.character(DF1$Genes)
    DF2$Genes<-as.character(DF2$Genes)

    DF2[DF2$Genes=='Unshared','Genes']<-paste(Prefix,'<',threshold[1],unit.DF1)
    levelOrderDF2<-factor(c('Identical','Similar','Mixed','Different',paste(Prefix,'<',threshold[1],unit.DF1)))
    if(length(threshold>1)){
      DF1[DF1$Genes=='Unshared','Genes']<-paste(Prefix,'<',threshold[2],unit.DF2)
      levelOrderDF1<-factor(c('Identical','Similar','Mixed','Different',paste(Prefix,'<',threshold[2],unit.DF2)))
    }else{
      DF1[DF1$Genes=='Unshared','Genes']<-paste(Prefix,'<',threshold[1],unit.DF2)
      levelOrderDF1<-factor(c('Identical','Similar','Mixed','Different',paste(Prefix,'<',threshold[1],unit.DF2)))
    }
    DF1$Genes<-factor(DF1$Genes,levelOrderDF1)
    DF2$Genes<-factor(DF2$Genes,levelOrderDF2)
  }
  if(P1){
    p1 <- ggplot(DF1,aes(x=factor(nb.tissues),colour=Genes,fill=Genes))
    p1 <- p1 + geom_bar(alpha=0.65)
    p1 <- p1 + labs(x='Number of tissues')

    if(annotate) p1 <- p1 + geom_text(stat='count',aes(label=..count..),position='stack',alpha=1,show.legend=FALSE)
    if(publish) p1 <- p1 + theme_bw(...)
    if(!missing(colour)) {
      p1 <- p1 + scale_colour_manual(values=colour)+scale_fill_manual(values=colour)
    }else{
      p1<- p1 + scale_color_brewer(palette = 'Set2')+scale_fill_brewer(palette='Set2')
    }
    print(p1)
  }
  if(P2){
    p2 <- ggplot(DF2,aes(x=factor(nb.tissues),colour=Genes,fill=Genes))
    p2 <- p2 + geom_bar(alpha=0.65)
    p2 <- p2 + labs(x='Number of tissues')

    if(annotate)  p2 <- p2 + geom_text(stat='count',aes(label=..count..),position='stack',alpha=1,show.legend=FALSE)
    if(publish) p2 <- p2 + theme_bw(...)

    if(!missing(colour)) {
      p2 <- p2 + scale_colour_manual(values=colour)+scale_fill_manual(values=colour)
    }else{
      p2 <- p2 + scale_color_brewer(palette = 'Set2')+scale_fill_brewer(palette='Set2')
    }
    print(p2)
  }
  #return(list(p1,p2))
  if(P1){
    if(P2){
      return(list(DF1,DF2,p1,p2))
    }else{
      return(list(DF1,DF2,p1))
    }
  }else{
    if(P2){
      return(list(DF1,DF2,p2))
    }else{
      return(list(DF1,DF2))
    }
  }
}


#' Compute and compare the breadth of expression of two datasets and send back the resuls as a list of data.frames
#'
#' @param DF1 numeric data.frame comprising the expression data of the first study.
#' @param DF2 numeric data.frame comprising the expression data of the second study.
#' @param omit.zero boolean. Default: TRUE.
#'                  Whether the genes which do not reach the minimal expression threshold in any tissue
#'                  should be stripped from the output.
#' @param threshold numeric. Default: 0. Minimal level of expression to be considered as expressed.
#' @param strict boolean. Default:FALSE. Whether if only strict equality between the breadth of expression of each gene
#'               between the two studies should be considered
#'               or if the "lapse" argument should be used to define 'similar' cases as well.
#' @param simplify boolean. Default: TRUE. Whether to focus only on the most extreme breadth of expression.
#'                 "Expressed in all" or "nearly in all" might be more descriptive than "expressed in 5" versus "expressed in 6".
#' @param Lims  a length-4 numeric vector. Default: c(1,3,10,12).
#'             When the analysis is in the "simplify" mode,
#'             it allows to delimit two ranges of expression breadth to compare between the two studies.
#' @param lapse positive integer. Default: 3. Allow to relax constraints on perfect equality of breadth for the genes between the two studies.
#'              A "similar" class is created. The lapse allows to define which level of similarity is acceptable.
#'              e.g. if lapse=3, all genes that have a breadth of expression that varies at most of 3 are considered
#'              to have similar breadth.
#'
#' @return a list of data.frames
#' @export
#'
xsharedDistrib<-function(DF1,DF2,omit.zero=TRUE,threshold=0,strict=FALSE,
                         simplify=TRUE,Lims=c(1,3,10,12),lapse=3){

  if(length(threshold)==1){
    DF1<-computeBreadth(DF1,threshold=threshold,omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold,omit.zero=omit.zero,typeR='df')
  }else{
    DF1<-computeBreadth(DF1,threshold=threshold[1],omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold[2],omit.zero=omit.zero,typeR='df')
  }

  if (missing(lapse)) lapse=3
  if(!simplify){
    DF1$Genes<-sharedBreadth(DF1,DF2,simplify=simplify,strict=strict,lapse=lapse)
    DF2$Genes<-sharedBreadth(DF2,DF1,simplify=simplify,strict=strict,lapse=lapse)
  }else{
    if(missing(Lims)) Lims=c(1,3,10,12)
    DF1$Genes<-sharedBreadth(DF1,DF2,simplify=simplify,strict=strict,Lims=Lims,lapse=lapse)
    DF2$Genes<-sharedBreadth(DF2,DF1,simplify=simplify,strict=strict,Lims=Lims,lapse=lapse)
  }

  return(list(DF1,DF2))
}


#' Plot the expression breadth distribution computed with xsharedDistrib
#'
#' @param DF a numeric data.frame that comprises the expression breadth of the genes
#' in a column named 'nb.tissues'
#' @param colour colour palette for the different categories.
#'               Default: "set2" of colorbrewer2
#' @param annotate boolean. Default: TRUE. Whether to annotate the figure
#'                 with the counts in each category for an accurate read.
#' @param publish  boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param ... other arguments that can be used by ggplot2::theme_bw()
#'
#' @return a figure
#' @export
#'
plot_xsharedDistrib<-function(DF,colour,annotate=TRUE,publish=TRUE,...){
  p <- ggplot(DF,aes(x=factor(nb.tissues),colour=Genes,fill=Genes))
  p <- p + geom_bar(alpha=0.5)
  p <- p + labs(x='Number of tissues')

  if(annotate) p <- p + geom_text(stat='count',aes(label=..count..),
                                  position='stack',alpha=1,show.legend=FALSE)

  if(publish) p <- p + theme_bw(...)

  if(!missing(colour)){
    p <- p + scale_colour_manual(values=colour)+scale_fill_manual(values=colour)
  }else{
    p <- p + scale_color_brewer(palette = 'Set2')+scale_fill_brewer(palette='Set2')
  }

  return(p)
}

#' Draw the figures of the shared expression for the two datasets.
#'
#' @param DF1 numeric data.frame comprising the expression data of the first study.
#' @param DF2 numeric data.frame comprising the expression data of the second study.
#' @param omit.zero boolean. Default: TRUE.
#'                  Whether the genes which do not reach the minimal expression threshold in any tissue
#'                  should be stripped from the output.
#' @param threshold numeric. Default: 0. Minimal level of expression to be considered as expressed.
#' @param strict boolean. Default:FALSE. Whether if only strict equality between the breadth of expression of each gene
#'               between the two studies should be considered
#'               or if the "lapse" argument should be used to define 'similar' cases as well.
#' @param colour vector of character strings to personalise the colours of the plot.
#' @param annotate boolean. Default: TRUE. Whether to annotate the figures with the counts in each category for an accurate read.
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param indexCol vector of integer. Comprise the index of the columns for which the expression breadth has to be compared.
#' @param ... more arguments to be treated
#'
#' @return a list of two figures.
#' @export
#'
cross.sharedDistrib_firstLast<-function(...,DF1,DF2,omit.zero=TRUE,threshold=0,strict=FALSE,
                                        colour,annotate=TRUE,publish=TRUE,indexCol){

  if(length(threshold)==1){
    DF1<-computeBreadth(DF1,threshold=threshold,omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold,omit.zero=omit.zero,typeR='df')
  }else{
    DF1<-computeBreadth(DF1,threshold=threshold[1],omit.zero=omit.zero,typeR='df')
    DF2<-computeBreadth(DF2,threshold=threshold[2],omit.zero=omit.zero,typeR='df')
  }
  reassign<-FALSE
  if(missing(indexCol)){
    indexCol=c(1,2,ncol(DF1))
    reassign<-TRUE
  }
  DF1$Genes<-sharedBreadth_firstLast(a=DF1,b=DF2,indexCol=indexCol,...)
  if(reassign) indexCol=c(1,2,ncol(DF2))
  DF2$Genes<-sharedBreadth_firstLast(a=DF2,b=DF1,indexCol=indexCol,...)


  p1 <- ggplot(DF1,aes(x=factor(nb.tissues),colour=Genes,fill=Genes))+geom_bar(position='stack',alpha=0.5)
  p1 <- p1 + labs(x='Number of tissues')
  if(annotate)
    p1 <- p1 + geom_text(stat='count',aes(label=..count..,
                                          y=ifelse(..count..>5,..count..,..count..*0.5)
    ),alpha=1,show.legend=FALSE)
  if(publish) p1 <- p1 + theme_bw()

  if(!missing(colour)) p1 <- p1 + scale_colour_manual(values=colour)+scale_fill_manual(values=colour)
  print(p1)

  p2 <- ggplot(DF2,aes(x=factor(nb.tissues),colour=Genes,fill=Genes,alpha=0.7))+geom_bar(position='stack',alpha=0.5)
  p2 <- p2 + labs(x='Number of tissues')
  if(annotate)
    p2 <- p2 + geom_text(stat='count',aes(label=..count..),alpha=1,show.legend=FALSE)

  if(publish) p2 <- p2 + theme_bw()

  if(!missing(colour)) p2 <- p2 + scale_colour_manual(values=colour)+scale_fill_manual(values=colour)
  print(p2)

  return(list(p1,p2))
}



#' Plot the genes in common between 3 datasets
#'
#' @param a numeric data.frame comprising the first dataset expression data
#' @param b numeric data.frame comprising the second dataset expression data
#' @param c numeric data.frame comprising the third dataset expression data
#' @param name a length-3 character string vector. Names of the three datasets.
#' @param print boolean. Default: TRUE. Whether to print the figure or not.
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param out character string. Default: "plot".
#'            Allows to pick the object returned by the function.
#' @param colorpalette vector of character strings to personalise the colours of the plot.
#' @param thm function of the type of theme(...)
#' @param xlab character string. Allows to change the label of the x-axis.
#' @param ylab character string. Allows to change the label of the y-axis.
#' @param baseSize size of the font. Default: 12. Applied only if "publish"
#' @param baseFamily name of the font.
#'                   Default is "Linux Libertine". Applied only if "publish"
#'
#' @return depending on the type of "out" can be a data.frame,
#'         a plot or the data to create the plot
#' @export
#'
barplotCommonGenesUnique3D<-function(a,b,c,name,print=TRUE,
                                     publish=TRUE,out='plot',
                                     colorpalette,thm,xlab,ylab,baseSize=12,
                                     baseFamily="Linux Libertine"){

  commonCol<-Intersect(colnames(a),colnames(b),colnames(c))
  new<-Reduce(rbind,lapply(commonCol, function(x){
    A<-rownames(a[a[,x]>0,])
    B<-rownames(b[b[,x]>0,])
    C<-rownames(c[c[,x]>0,])
    abc <- Intersect(A,B,C)
    ab  <- setdiff(Intersect(A,B),abc)
    ac  <- setdiff(Intersect(A,C),abc)
    bc  <- setdiff(Intersect(B,C),abc)
    a   <- setdiff(A,Union(abc,ab,ac))
    b   <- setdiff(B,Union(abc,ab,bc))
    c   <- setdiff(C,Union(abc,ac,bc))
    rbind(data.frame(Tissue=x,ID=a,Group=paste(name[1],'only')),
          data.frame(Tissue=x,ID=b,Group=paste(name[2],'only')),
          data.frame(Tissue=x,ID=c,Group=paste(name[3],'only')),
          data.frame(Tissue=x,ID=ab,Group=paste(name[1],"&",name[2])),
          data.frame(Tissue=x,ID=ac,Group=paste(name[1],"&",name[3])),
          data.frame(Tissue=x,ID=bc,Group=paste(name[2],"&",name[3])),
          data.frame(Tissue=x,ID=abc,Group='All 3 datasets')

    )
  }))
  p <- ggplot(new, aes(x=Tissue))+geom_bar(aes(y=..count..,fill=Group))+coord_flip()
  if(publish) p<- p+theme_bw(base_size=baseSize,base_family = baseFamily)
  if(!missing(thm)) p<-p+thm
  #p <- p + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
  p <- p + theme(legend.position = 'bottom')
  p <- p + guides(fill=guide_legend(title="Present in",reverse=TRUE))
  if(!missing(colorpalette)){
    #  if(length(colorpalette)<length(unique(new$Group))){
    #    names(colorpalette)<-unique(new$Group)
    p <- p + scale_fill_manual(values=colorpalette,drop=FALSE)
    #  }else{
    #    message("non adequate colorpalette (size)")
    #  }
  }
  if(!missing(xlab))  p<-p+xlab(xlab)
  if(!missing(ylab))  p<-p+labs(y=ylab)
  if(print) print(p)

  switch(out,
         'DF'            = return(new),
         'plot'          = return(print(p)),
         'plot_data'     = return(p))

}


#' Plot the genes in common between 2 datasets
#'
#' @param a numeric data.frame comprising the first dataset expression data
#' @param b numeric data.frame comprising the second dataset expression data
#' @param name a length-2 character string vector. Names of the three datasets.
#' @param print boolean. Default: TRUE. Whether to print the figure or not.
#' @param publish boolean. Default: TRUE. Whether to apply ggplot2::theme_bw to the plot.
#' @param out character string. Default: "plot".
#'            Allows to pick the object returned by the function.
#' @param colorpalette vector of character strings to personalise the colours of the plot.
#' @param thm function of the type of theme(...)
#' @param xlab character string. Allows to change the label of the x-axis.
#' @param ylab character string. Allows to change the label of the y-axis.
#' @param baseSize size of the font. Default: 12. Applied only if "publish"
#' @param baseFamily name of the font.
#'                   Default is "Linux Libertine". Applied only if "publish"
#'
#' @return depending on the type of "out" can be a data.frame,
#'         a plot or the data to create the plot
#' @export
#'
barplotCommonGenesUnique2D<-function(a,b,name,print=TRUE,
                                     publish=TRUE,out='plot',
                                     colorpalette,thm,xlab,ylab,baseSize=12,
                                     baseFamily="Linux Libertine"){

  commonCol<-base::intersect(colnames(a),colnames(b))
  new<-Reduce(rbind,lapply(commonCol, function(x){
    A<-rownames(a[a[,x]>0,])
    B<-rownames(b[b[,x]>0,])
    ab <- Intersect(A,B)
    ua <- setdiff(A,ab)
    ub <- setdiff(B,ab)

    rbind(data.frame(Tissue=x,ID=ua,Group=paste(name[1],'only')),
          data.frame(Tissue=x,ID=ub,Group=paste(name[2],'only')),
          data.frame(Tissue=x,ID=ab,Group=paste('Both',name[1],'&',name[2]))
    )
  }))
  p <- ggplot(new, aes(x=Tissue))+geom_bar(aes(y=..count..,fill=Group))+coord_flip()
  if(publish) p<- p+theme_bw(base_size=baseSize,base_family = baseFamily)
  if(!missing(thm)) p<-p+thm
  p <- p + theme(legend.position = 'bottom')
  p <- p + guides(fill=guide_legend(title="Present in",reverse=TRUE))
  if(!missing(colorpalette)){
    p <- p + scale_fill_manual(values=colorpalette,drop=FALSE)
  }
  if(!missing(xlab))  p<-p+xlab(xlab)
  if(!missing(ylab))  p<-p+labs(y=ylab)
  if(print) print(p)

  switch(out,
         'DF'            = return(new),
         'plot'          = return(print(p)),
         'plot_data'     = return(p))

}


#' Returns the id (rownames) of the genes that are found in the tissues that have a
#' breadth of expression equals to or lesser than a given one
#'
#' @param DFa data.frame; should be logical can be numeric if comp is TRUE
#' @param nbMax integer; breadth of expression maximal to consider
#' @param naCol character string. Name of the column which gives the breadth of expression of each gene
#'              default: "nb.tissues"
#' @param comp logical; should the breadth of expression be computed first from a numeric DFa
#' @param threshold numeric, which cutoff of expression has to be used for considering a gene expressed in a tissue or not
#' @param omit.zero logical; whether the rows with zeros only should be kept or not
#'
#' @return a list that comprises for each column of the inputted data.frame,
#'         a vector with the rownames of the genes expressed in a maximal number of tissues (given).
#' @export
#'
ExpressedGeneInTissue<-function(DFa,nbMax=2,
                                naCol="nb.tissues",
                                comp=FALSE,threshold,omit.zero){

  if(comp) DFa<-computeBreadth(DFa,
                                  threshold=threshold,
                                  omit.zero=omit.zero,typeR='dfBool')

  DFa<-DFa[DFa[,naCol]<nbMax+1,]
  DFa<-DFa[,setdiff(colnames(DFa),naCol)]

  DFList<-lapply(setNames(colnames(DFa),colnames(DFa)),
                    function(x){
                      return(rownames(DFa)[DFa[,x]])
                    })

  return(DFList)
}


