# Retrieve external data, e.g. from biomaRt --------------

#' Retrieve gene symbols (hgnc) for Ensembl gene IDs with biomaRt
#' @description function which has generated "gene.mapID"
#'
#' @param geneID vector of Ensembl gene ids
#' @param out 'all' for a named vector where the gene symbol are the content
#'                  and the Ensembl ids are the names
#'             'DF' for a data.frame where the first column contains the Ensembl ids
#'             and the second one the gene symbols.
#'             'nameOnly' for the gene symbols only
#' @param ens Formal class 'Mart' from package biomaRt
#'
#' @return gene symbols (with or without their corresponding Ensembl ids)
#' @export
#'
nameID<-function(geneID,out='nameOnly',ens){
  if(missing(ens)){
    ens = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  }

  new<-biomaRt::getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                      mart=ens,filters="ensembl_gene_id",values=geneID)
  switch(out,
         'all' = {return(setNames(new$hgnc_symbol,new$ensembl_gene_id))},
         'DF'  = {return(new)},
         'nameOnly'= return(new[new$ensembl_gene_id %in%geneID,]$hgnc_symbol)
  )
}


#' Retrieve additional gene information from biomaRt
#'
#' @param geneID Ensembl gene ids to query
#' @param ens a "Mart" object from biomaRt
#' @param output "DF" for dataframe,"data.table" for a data.table::data.table,widget for a DT::datatable
#'               and "latex" for a latex formatted output
#' @param otherdata other data that can be queried with biomaRt::getBM.
#'                  See possible attributes=biomaRt::listAttributes(ens)
#' @param NoIDcol logical; default=TRUE: remove the first column that contains the ensembl gene id
#'
#' @return a data.frame or a data.table or latex code or widget
#' @export
#'
tableGeneIDdesc<-function(geneID,ens,output='DF',otherdata=NULL,NoIDcol=TRUE){
  if(missing(ens)){
    ens = biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")
  }

  queryAtt<-c("ensembl_gene_id",
              "hgnc_symbol",
              "description",
              "gene_biotype")

  if(!missing(otherdata)) queryAtt<-c(queryAtt,otherdata)

  new<-biomaRt::getBM(attributes=queryAtt,
                      mart=ens,filters="ensembl_gene_id",
                      values=geneID)

  rownames(new)<-new$ensembl_gene_id
  colnames(new)<-c('Ensembl gene ID','Gene name','Description', 'Gene biotype')

  if(NoIDcol) new<-new[,-1]

  switch(output,
         "widget"     = {DT::datatable(new)},
         "DF"         = {return(new)},
         "data.table" = {return(data.table::data.table(new))},
         "latex"      = {return(xtable::xtable(new,booktabs=TRUE))}
  )
}


#' Retrieve GO information for a list of genes or all genes
#'
#' @param queryList character vector; optional; id of genes
#' @param biomart default: 'ENSEMBL_MART_ENSEMBL'
#' @param dataset default: 'hsapiens_gene_ensembl'
#' @param attributes default: c("ensembl_gene_id","name_1006", "go_id")
#' @param filters default; "ensembl_gene_id"
#' @param host optional;
#'
#' @return a data.frame with go identifiers and their names.
#' @export
#'
createGoAnnot<-function(queryList,
                        biomart='ENSEMBL_MART_ENSEMBL',
                        dataset='hsapiens_gene_ensembl',
                        attributes=c("ensembl_gene_id",
                                     "name_1006",
                                     "go_id"),
                        filters=("ensembl_gene_id"),
                        host){

  if(!missing(host)){
    ensembl<-biomaRt::useMart(host=host, biomart=biomart, dataset=dataset)
  }else{
    ensembl<-biomaRt::useMart(biomart=biomart, dataset=dataset)
  }

  if(missing(queryList)){
    annotGO<-biomaRt::getBM(attributes=attributes, mart=ensembl)
  }else{
    annotGO<-biomaRt::getBM(attributes=attributes,mart=ensembl,
                            filters=filters,values=queryList)
  }
  return(annotGO)
}

#' Takes a list of GO ids and queries all the genes related to it (including the ones for child nodes)
#'
#' @param goID list of GO ids to retrieve the genes for
#' @param annotGO optional; data.frame with information of the go ids and the related genes
#' @param pathGeneMap path to the file for the annotGO
#' @param queryList argument for createAnnot if no annotGO or pathGeneMap are not provided
#' @param GOList GO ids related to the list of the provided GO ids (such as their child nodes)
#' @param ontology one from 'BP','CC' or 'MF'
#' @param exportAnnotGO logical; whether annotGO should be exported to the environement sharedenv for future use
#' @param exportGOList logical; whether GOList should be exported to the environement sharedenv for future use
#' @param verbose logical; allows to better track how annotGO and GOList are created/used
#' @param ... other argument to pass to createGoAnnot
#'
#' @return Gene list associated to the GO ids
#' @export
#'
go2genes<-function(goID,annotGO,
                   pathGeneMap,
                   queryList,
                   GOList,ontology,
                   exportAnnotGO=TRUE,
                   exportGOList=TRUE,
                   verbose=FALSE,
                   ...){

  if(missing(goID)) return("goID is a non-optional argument")
  if(missing(annotGO)&&missing(pathGeneMap)){
    annotGO<-try(get('annotGO',envir=sharedEnv),silent=TRUE)
    if(class(annotGO)=="try-error") {
      if(!missing(queryList)){
        annotGO<-try(createGoAnnot(queryList=queryList,...),silent=TRUE)
      }
      annotGO<-try(createGoAnnot(...),silent=TRUE)
      if(class(annotGO)=="try-error") {
        utils::data("goAnnot.76")
        message("Warning: createGoAnnot failed; used data('goAnnot.76') for annotGO")
      }
    }
  }else{
    if(!missing(pathGeneMap)){
      if(file.exists(pathGeneMap)){
        annotGO<-try(Load(pathGeneMap),silent=TRUE)
        if(class(annotGO)=="try-error") {
          utils::data("goAnnot.76")
          message(paste("Warning: loading content of",
                        pathGeneMap,"failed; used data('goAnnot.76') for annotGO"))
        }
      }else{
        message('you need to provide a valid file path')
      }
    }
  }
  if(exportAnnotGO) assign('annotGO',value=annotGO,envir=sharedEnv)

  if(missing(GOList)){
    GOList<-try(get('GOList',envir=sharedEnv),silent=TRUE)
    if(class(GOList)=="try-error") {
      if(!exists("GOSimEnv")){
        if(!missing(ontology)&& ontology %in% c("BP","MF","CC")){
          GOSim::setOntology(ont = ontology)
        }else{
          GOSim::setOntology()
        }
      }else{
        if(!missing(ontology)&&!ontology==GOSimEnv$ontology){
          GOSim::setOntology(ont = ontology)
        }
      }
      GOList<-GOSim::getOffsprings()
    }
  }

  #For each GO term given, the list of GO terms childs (even indirect ones is retrieved)
  #then from the mapping between the GO terms and the genes all the genes that have a GO term
  # that is from the list we want (even if this is an indirect go term) are retrieved
  geneList<-c(unlist(lapply(goID,function(x){
    lapply(GOList[[x]],function(y){
      annotGO[annotGO$go_id==y,]$ensembl_gene_id})
  })
  ))

  if(exportGOList) assign('GOList',value=GOList,envir=sharedEnv)

  return(unique(geneList))

}
utils::globalVariables("GOSimEnv")

# Plots and other outputs --------------------------------------------

#' Create a heatmap from a compareClusterResult (\code{\link[clusterProfiler]{compareCluster}}) package)
#' showing for each enriched GO term, the overall observed correlation
#' of the associated genes set.
#'
#' @param compareCluster Result object of \code{\link[clusterProfiler]{compareCluster}}
#' @param level integer; default: 5. How many level of the comparison between cluster should be retrieved
#' @param corrValue numeric vector. The correlation values for the different GO terms associated gene sets.
#' @param corr logical; default: TRUE. Whether the correlation should be plotted.
#' @param sortDF logical; default: TRUE. Whether the results should be sorted by correlation.
#' @param simpleHeatmap logical; default: TRUE.
#'                      Whether to remove the GO for which the gene set is empty.
#' @param corName string; default: "Pearson Correlation".
#'                Allow to personalise the correlation legend
#' @param base_family string; default: "Linux Libertine". Name of the font family to use for the plot
#' @param base_size numeric; size of the main text of the plot; default: 16
#' @param stripName logical; default: FALSE.
#' @param out string; switch that allows to chose the type result object.
#'            default:  "DF", gives a data.frame (used for the plot),
#'            other options: "mat" gives the previous result formatted as a matrix
#'                            "plot" gives the figure ready to be saved as an object.
#' @param plot logical; default: TRUE. Whether a plot should outputed to the standard output.
#'
#' @return depends of the value of "out"
#' @export
#'
plotGOHeatmap<-function(compareCluster,level=5,
                        corrValue,corr=TRUE,sortDF=TRUE,
                        simpleHeatmap=TRUE,
                        corName="Pearson\nCorrelation",
                        base_family='Linux Libertine',
                        base_size=16,stripName=FALSE, out='DF',
                        plot=TRUE){

  extract<-clusterProfiler::compareCluster@compareClusterResult

  goid<-parallel::mclapply(levels(extract$Cluster),function(x){
    if(length(extract[extract$Cluster== x,'ID'])>level){
      return(extract[extract$Cluster== x,'ID'][1:level])
    }else{
      return(extract[extract$Cluster== x,'ID'])
    }})
  DFres<-parallel::mclapply(goid,function(x){
    tmp2<-parallel::mclapply(x,function(y){
      tags<-unlist(strsplit(extract[extract$ID==y,]$geneID,'/'))
      if(exists('tags')){
        tempo<-as.integer(names(corrValue) %in% tags)
      }else{
        tempo<-corrValue*0
      }
      rm(tags)
      if(corr) tempo<-tempo*corrValue
      return(tempo)})
    return(tmp2)})


  DFres<-data.frame(DFres)
  colnames(DFres)<-unlist(goid)
  rownames(DFres)<-names(corrValue)

  DFtoReturn<-DFres
  DFres<-data.frame('Correlation'=corrValue,DFres)

  if(simpleHeatmap) { #keep only genes and GO terms with no null value
    DFres<-DFres[rowSums(DFres[,-1]!=0)>0,]
    DFres<-DFres[,colSums(DFres)!=0]
  }

  DFres<-DFres[ rownames(DFres[order(DFres[,1],decreasing=TRUE),]) ,]

  colnames(DFres)<-c('Correlation',
                     extract[gsub(':','.',extract$ID) %in% colnames(DFres),"Description"])

  melted_t_DFres<-melt(t(DFres))
  melted_t_DFres$Var1<-factor(melted_t_DFres$Var1,
                              levels = rev(levels(melted_t_DFres$Var1)))

  p <- ggplot(data = melted_t_DFres,
              aes_string(x="Var2", y="Var1", fill="value")) + geom_tile(color = "white")
  p <- p + theme_minimaliste(base_family=base_family,base_size=base_size)
  p <- p + scale_fill_gradient2(low = "blue", high = "red", mid = "white",
                                midpoint = 0, limit = c(-1,1), space = "Lab",
                                name=corName)
  p <- p + theme(legend.position = 'left')+xlab("")+ylab("")
  p <- p + theme(axis.title.x=element_blank(),
                 axis.text.x=element_blank(),
                 axis.ticks.x=element_blank())
  p <- p + theme(plot.margin=unit(c(0,0,-1,0), "cm"),
                 legend.margin=margin(0,0,0,0),
                 legend.box.margin=margin(0,0,0,0))

  if(stripName)  p <- p+theme(axis.title.y=element_blank(),
                              axis.text.y=element_blank(),
                              axis.ticks.y=element_blank())

  if(plot) print(p)

  if(sortDF) DFtoReturn<-DFtoReturn[ rownames(DFtoReturn[order(DFtoReturn[,1],
                                                               decreasing=TRUE),]) ,]
  switch (out,
          "DF" =  return(DFtoReturn),
          "mat"=  return(as.matrix(DFtoReturn)),
          "plot"= return(p)
  )
}

