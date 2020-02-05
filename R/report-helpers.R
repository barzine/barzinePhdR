# layout helpers------------------

#' Create a customised HTML table widget based on the DT library
#'
#' @param data a data object (either a matrix or a data frame)
#' @param ... other parameters handled by DT::datatable
#'
#' @return an HTML table widget
#' @export
#'
datatable.custom<-function(data,...){
  DT::datatable(data, ...,
                extensions= c('Buttons','ColReorder'),
                options = list(dom='l<"clear">f<"clear">Brtip',colReorder = TRUE,
                               buttons = list('copy','print',
                                              list(extend='collection',
                                                   buttons=c('csv', 'excel', 'pdf'),
                                                   text='Download'),
                                              I('colvis')),
                               pageLength = 15,
                               pageLength = 15,
                               lengthMenu = c(5, 10, 20,50,100)
                ))
}


#' Wrapper around knitr::kable with HTML as a default
#'
#' @param ... any parameter that can be treated by knitr::kable
#' @param format default format="HTML"
#'
#' @return a character vector of the table source code
#' @export
#'
#' @seealso \code{\link[knitr]{kable}}
#'
cleantable<-function(...,format="html"){
  knitr::kable(...,format=format)
}

#' Encapsulate result between <pre></pre> html tags
#'
#' @param string content to be printed between the <pre></pre> tags
#'
#' @return a bloc of characters
#' @export
#'
cat.html<-function(string){
  cat("\n\n<pre>\n")
  print(string)
  cat("\n</pre>\n\n")
}



#' Format results output for two-dataframe meta-analysis
#'
#' @param DF1 first data.frame to include
#' @param DF2 second data.frame to include
#' @param GeneNameDic named character vector that maps the gene names to their identifiers
#' @param externalSupport allows to include other information about the genes
#' @param select 'all_common' to keep only the genes that are found in both data.frames
#'               'all' to include all genes from both data.frames
#' @param MEAN logical; default: TRUE. Gives the mean value of the gene expression
#' @param SD logical; default: TRUE. Gives the standard deviation of the expression of the genes.
#' @param CV logical; default: TRUE. Gives the coefficient variation of the expression of the genes
#' @param nb.Tissues logical; default: TRUE. The expression breadth of the genes
#'                   (i.e. how many tissues the genes has been detected)
#' @param MEDIAN logical; default: TRUE. Median expression of the genes
#' @param MAD logical; default: TRUE median absolute deviation of the expression of the genes
#' @return a data.frame summarising different aspect of the expression of the genes
#' @export
#'
createTableResults_basics<-function(DF1,DF2,GeneNameDic,externalSupport,select='all_common',
                                    MEAN=TRUE,SD=TRUE,CV=TRUE,nb.Tissues=TRUE,
                                    MEDIAN=TRUE,MAD=TRUE){

  nameDF1<-deparse(substitute(DF1))
  nameDF2<-deparse(substitute(DF2))

  switch(select,
         'all'={
           geneList<-union(rownames(DF1),rownames(DF2))
         },
         'all_common'={
           geneList=intersect(rownames(DF1),rownames(DF2))
         },
         {
           geneList=intersect(rownames(DF1),rownames(DF2))
         }
  )

  DF1<-DF1[geneList,]
  DF2<-DF2[geneList,]

  ## Processing data
  toBind<-list(Gene.name=GeneNameDic[geneList],
               Pearson.corr=compute_gene_cor(DF1,DF2,'pearson'),
               Pearson.corr.Log2=compute_gene_cor(DF1,DF2,'pearson_log2'),
               Spearman.corr=compute_gene_cor(DF1,DF2,'spearman'))

  if(nb.Tissues){
    toBind[[paste0('nb.Tissues.',nameDF1)]]<-rowSums(DF1>0)
    toBind[[paste0('nb.Tissues.',nameDF2)]]<-rowSums(DF2>0)
  }

  if(MEAN){
    toBind[[paste0('Mean.',nameDF1)]]<-rowMeans(DF1)
    toBind[[paste0('Mean.',nameDF2)]]<-rowMeans(DF2)
  }

  if(SD){
    toBind[[paste0('sd.',nameDF1)]]<-rowSD(DF1)
    toBind[[paste0('sd.',nameDF2)]]<-rowSD(DF2)
  }

  if(CV){
    toBind[[paste0('cv.',nameDF1)]]<-rowSD(DF1)/rowMeans(DF1)
    toBind[[paste0('cv.',nameDF2)]]<-rowSD(DF2)/rowMeans(DF2)
  }

  if(MEDIAN){
    toBind[[paste0('Median.',nameDF1)]]<-RowMedians(DF1)
    toBind[[paste0('Median.',nameDF2)]]<-RowMedians(DF2)
  }

  if(MAD){
    toBind[[paste0('mad.',nameDF1)]]<-RowMad(DF1)
    toBind[[paste0('mad.',nameDF2)]]<-RowMad(DF2)
  }

  if(!missing(externalSupport)) {
    if(is.vector(as.character(externalSupport))){
      toBind[['Evidence']]<-as.character(externalSupport[geneList])
    }
  }

  keltype<-sapply(toBind,class)
  DFres<-do.call(cbind,toBind)
  DFres<-as.data.frame(DFres,stringsAsFactors=FALSE)
  ColDFres<-colnames(DFres)
  DFres<-data.frame(lapply(colnames(DFres),function(x){
    return(do.call(paste0('as.',keltype[[x]]),list(DFres[,x])))
  }),
  stringsAsFactors=FALSE)
  colnames(DFres)<-ColDFres
  rownames(DFres)<-geneList
  return(DFres)
}

#' Wrapper automating the creation of scatter plots (scatplot.log2)
#' based on scatplot.log2
#'
#' @param DFindex When used for the output for "identical" (t.test output) should be already sorted
#' @param DF1 first numeric data.frame
#' @param DF2 second numeric data.frame
#' @param label.df1 name on the title for DF1 (and used to create the filename)
#' @param label.df2 name on the title for DF2 (and used to create the filename)
#' @param out.df1 name to output on the x-axis (for DF1)
#' @param out.df2 name to output on the y-axis (for DF2)
#' @param xmin lower limit for the x-axis
#' @param xmax upper limit for the x-axis
#' @param ymin lower limit for the y-axis
#' @param ymax upper limit for the y-axis
#' @param a numeric; alpha scale
#' @param xcor numeric; x coordinate to center the correlation
#' @param ycor numeric; y coordinate to center the correlation
#' @param methodcor "spearman" or "pearson"
#' @param fig switch; "grob"; "all"
#' @param abline logical; default: FALSE. Whether x=y should be drawn.
#' @param verbose logical; default: TRUE.
#' @param publi logical; default: TRUE
#' @param report logical; default: TRUE
#' @param pdf logical; default: FALSE
#'
#' @return several scatter plots outputted as a pdf or part as html report
#' @export
#'
wrap.scatplot.log2<-function(DFindex,DF1,DF2,label.df1,label.df2,out.df1,out.df2,
                             xmin,xmax,ymin,ymax,a,xcor,ycor,methodcor,fig,
                             abline=FALSE,verbose=TRUE,publi=TRUE,
                             report=TRUE,pdf=FALSE){

  if (report) cat(paste0("\n\n### Scatter plots for ",label.df1," and ", label.df2,"\n\n"))
  colnames(DFindex)[1:2]<-c('DF1','DF2')

  if(pdf) fig='grob'

  if(missing(fig)) fig='all'

  tmp<-parallel::mclapply(1:nrow(DFindex),function(numb){
    if(report) cat(paste0("\n\n#### ",label.df2,": ",as.character(DFindex$DF2[numb]),
                          " ~ ",label.df1,": ",as.character(DFindex$DF1[numb]),"\n\n"))
    DF<-data.frame(DF1=DF1[,as.character(DFindex[numb,'DF1'])],
                   DF2=DF2[,as.character(DFindex[numb,'DF2'])])
    rownames(DF)<-rownames(DF1)
    p.tmp<-scatplot.log2(DF,'DF2','DF1',title=paste(label.df2,':',as.character(DFindex$DF2[numb]),
                                                    '~',label.df1,':',as.character(DFindex$DF1[numb])),
                         xlabs=paste(out.df2,':',as.character(DFindex$DF2[numb])),
                         ylabs=paste(out.df1,':',as.character(DFindex$DF1[numb])),
                         xmin=xmin, xmax=xmax,ymin=ymin,ymax=ymax,a=a,
                         xcor=xcor,ycor=ycor,methodcor=methodcor,
                         abline=abline,verbose=verbose,
                         publi=publi,fig=fig)
    if(!pdf){
      print(p.tmp)
    }else{
      ggsave(filename =
               paste0('scatterplot_',label.df1,'_',label.df2,'_',
                      as.character(DFindex$DF2[numb]),'.pdf'),
             p.tmp)
    }
    rm(DF)
    return()
  })
}


