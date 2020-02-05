## Generic functions -------------------

#' Change the name of the samples based on preset annotation tables
#' @description Specific to the phd projet
#'
#' @param DF a data.frame that the columns names have to be changed
#' @param originalNames the original column names
#' @param newNames the new column names
#' @param annot.DF if it exists, takes over originalNames and newNames
#' @param prefix string will be add to the colnames after the change of id
#' @param suffix string will be removed from the data.frame colnames before the change of id
#'
#' @return a data.frame with changed colnames
#' @export
#'
change.localID<-function(DF,originalNames,newNames,annot.DF, prefix, suffix){
  if(!missing(suffix)){
    coltmp<-sub(suffix,'',colnames(DF))
    colnames(DF)<-coltmp
  }
  if(!missing(originalNames)&&!missing(newNames)&&missing(annot.DF)){
    if (length(originalNames)==length(newNames)){
      sample<-c(originalNames)
      localID<-c(newNames)
      annot.DF<-data.frame(sample,localID)
    }else{
      if(length(originalNames)<length(newNames))
        print("ERROR: the current identifiers have too be as many as the new identifiers")
      if(length(originalNames)>length(newNames))
        print("ERROR: the new identifiers have too be as many as the current identifiers")
      return(print('operation failed'))
    }
  }

  if(!missing(annot.DF))
  {
    message("Only the information contains in the annotation data.frame is used")
  }

  if(!missing(prefix))
  {
    annot.DF$localID<-sapply(annot.DF$localID, function(x) paste(prefix,x,sep='.'))
  }

  colnamesTemp<-sapply(colnames(DF),function(x){
    as.character(annot.DF[annot.DF$sample==x,]$localID)
  }
  )
  colnames(DF)<-colnamesTemp
  return(DF)
}


#' Merge the replicates of a data.frame together
#'
#' @param annot.DF annotation dataframe that indicates
#'                 the organisation of the biological and technical replicates
#' @param DF data.frame
#' @param group "double" (default); "Biological" ("b" or "B"); "Technical" ("t" or "T")
#' @param FUN the function with which the replicates (biological and/or technical) are merged. Default: "rowMeans"
#' @param FUNT function with which the technical replicates are merged
#' @param FUNB function with which the biological replicates are merged
#' @param ... not implemented
#'
#' @return data.frame
#' @export
#'
merge_replicates<-function(annot.DF, DF,group="double",FUN='rowMeans',FUNT,FUNB,...){
  if (missing(annot.DF))
    stop("Annotation missing!")
  if (missing(DF))
    stop("Data missing!")

  switch(group,
         't'=,
         'T'=,
         'Technical'={
           if(missing(FUNT)){
             return(merge_tech_rep(annot.DF,DF,FUN,...))
           }else{
             return(merge_tech_rep(annot.DF,DF,FUNT,...))
           }
         },
         'b'=,
         'B'=,
         'Biological'={
           if(missing(FUNB)){
             return(merge_biol_rep(annot.DF,DF,FUN,...))
           }else{
             return(merge_biol_rep(annot.DF,DF,FUN,...))
           }

         },
         'double'={
           if(missing(FUNT)){
             new.DF<-merge_tech_rep(annot.DF,DF,FUN,...)
           }else{
             new.DF<-merge_tech_rep(annot.DF,DF,FUNT,...)
           }
           if(missing(FUNB)){
             return(merge_biol_rep(annot.DF,new.DF,FUN,...))
           }else{
             return(merge_biol_rep(annot.DF,new.DF,FUNB,...))
           }
         },
         {print('You need to define which replicates have to be pooled\n"T" or "Techninal" for technical replicates and "B" or "Biological" for biological replicates')
           return(NA)
         }
  )
}


#' Merge the techincal replicates in a data.frame together
#'
#' @param annotDF annotation dataframe that indicates the organisation of the technical replicates
#' @param DF data.frame
#' @param FUN the function with which the replicates technical are merged. Default: "rowMeans"
#'
#' @return a data.frame
#' @export
#'
merge_tech_rep<-function(annotDF,DF,FUN='rowMeans'){
  if (missing(annotDF))
    stop("Annotation missing!")
  if (missing(DF))
    stop("Data missing!")

  if(!"BioID"%in% colnames(annotDF)){
    print("Impossible to retrieve technical replicates: 'BioID' missing in annotation")
    return(DF)
  }

  new.DF<-as.data.frame(lapply(unique(annotDF$condition),function(x){
    ID<-unique(annotDF[annotDF$condition==x,"BioID"])
    tempL<-lapply(ID,function(y) {
      if(class(DF[,annotDF[annotDF$condition==x&annotDF$BioID==y,"sample"]])!="data.frame"){
        return(DF[,annotDF[annotDF$condition==x&annotDF$BioID==y,"sample"]])
      }else{
        return(call(FUN,DF[,annotDF[annotDF$condition==x&annotDF$BioID==y,"sample"]]))
      }

    })
    names(tempL)<-sapply(ID,function(sub){paste(tolower(x),sub,sep='.')})
    return(tempL)


  }))
  rownames(new.DF)<-rownames(DF)
  return(new.DF)
}


#' Merge the biological replicates in a data.frame together
#'
#' @param annotDF annotation dataframe that indicates the organisation of the biological replicates
#' @param DF data.frame
#' @param FUN the function with which the replicates biological are merged. Default: "rowMeans"
#' @param alone if all the columns are replicates of a single sample; default=FALSE
#' @param ... not implemented
#'
#' @return a data.frame
#' @export
#'
merge_biol_rep<-function(annotDF,DF,FUN='rowMeans',alone=FALSE,...){
  if(!alone){
    if (missing(annotDF))
      stop("Annotation missing!")
    if (missing(DF))
      stop("Data missing!")

    new.DF<-as.data.frame(lapply(unique(annotDF$condition),function(x){
      if(!is.data.frame(DF[,grep(x,colnames(DF),ignore.case=TRUE)])){
        return(DF[,grep(x,colnames(DF),ignore.case=TRUE)])
      }else{
        if(length(unique(annotDF[annotDF$condition==x,"lib.type"]))>1)
          warning("Several library types included in the merged output")
        if(length(unique(annotDF$condition[grep(x,annotDF$condition,ignore.case=TRUE)]))>1){
          allMatch<-colnames(DF)[grep(x,colnames(DF),ignore.case=TRUE)]
          otherMatchingPattern<-setdiff(annotDF$condition[grep(x,annotDF$condition,ignore.case=TRUE)],x)
          uncorectMatch<-colnames(DF)[unlist(lapply(otherMatchingPattern, function(p) {
            grep(p,colnames(DF),ignore.case=TRUE)}))]
          correctMatch<-setdiff(allMatch,uncorectMatch)
          return(call(FUN,DF[,correctMatch]))
        }else{
          # if(length(grep(x,colnames(DF),ignore.case=TRUE))>0){
          return(call(FUN,DF[,grep(x,colnames(DF),ignore.case=TRUE)]))
          # }#else{
          #      if(any(colnames(DF)%in%annotDF[annotDF$condition==x,]$sample)){
          #          return(do.call(FUN,list(as.data.frame(DF[,annotDF[annotDF$condition==x,]$sample]))))
          #      }
          #  }
        }
      }
    }
    ))
    colnames(new.DF)<-unique(annotDF$condition)
    rownames(new.DF)<-rownames(DF)
    return(new.DF)
  }else{
    new.DF<-data.frame(lapply(unique(annotDF$condition),function(x){
      if(class(DF[,annotDF[annotDF$condition==x,"sample"]])!="data.frame"){
        return(DF[,annotDF[annotDF$condition==x,"sample"]])
      }else{
        return(do.call(FUN,list(DF[,annotDF[annotDF$condition==x,]$sample])))
      }
    }))

    colnames(new.DF)<-unique(annotDF$condition)
    rownames(new.DF)<-rownames(DF)
    return(new.DF)
  }
}


#' Give the maximum value observed for each gene
#' across all samples of a given tissue/condition
#'
#' @param DF numeric data.frame
#' @param annot.df annotation data.frame
#' @param cond vector of conditions (tissues) to be found in the returned data.frame
#'
#' @return a data.frame
#' @export
#'
maxCond.df<-function(DF,annot.df,cond='all'){
  if (cond[1]=='all'){
    cond<-as.character(unique(annot.df$cond))
  }
  #population of the matrix
  listCond<-parallel::mclapply(cond, function(x) {
    cols<-as.character(annot.df[annot.df$condition==x,"sample"])
    if(length(cols)>1){
      apply(DF[,cols],1,max)
    }else{
      if(length(cols)==1){
        identity(DF[,cols])}   #else{
      #might need something
      #}
    }
  })
  df.max<-Reduce(cbind,listCond)
  colnames(df.max)<-cond
  rownames(df.max)<-rownames(DF)
  return(as.data.frame(df.max))
}

## More specific preliminaries for individual datasets -----------------------

#' Remove the proteins that are part of a cluster
#'
#' @param DF numeric (expression) data.frame
#' @param path path to the file that provide the clusters description
#'
#' @return a numeric data.frame where the problematic proteins have been removed
#' @export
#'
removeClusters<-function(DF,path){
  clusterList<-Read.table(path,row.names=NULL)[,1]
  clusterList<-clusterList[grep(',',clusterList)]
  clusterList<-unique(Reduce(c,strsplit(clusterList,',')))
  return(DF[!rownames(DF) %in% clusterList,])
}

#' Allows to create a function that extracts the columns
#' that contains a specific tag
#'
#' @param tag character string that should be look for
#'
#' @return a function
#' @export
#'
extractWithTag<-function(tag){
  function(DF,clean=TRUE){
    DF<-DF[,grep(tag,colnames(DF))]
    if(clean)  colnames(DF)<-gsub(tag,'',colnames(DF))
    return(DF)
  }
}

#' Extract columns that contains "WithInSampleTop3Abundance." in their name
#'
#' @param DF an expression data.frame
#' @param clean boolean; default: TRUE.
#'              Whether rows that contains only null value should be removed
#'
#' @return  a data.frame which contains only the columns that contains the tag
#' @export
#'
extractWithinSampleTop3Abundance<-extractWithTag(tag="WithInSampleTop3Abundance.")


#' Extract columns that contains "WithInSampleAbundance." in their names
#'
#' @param DF an expression data.frame
#' @param clean boolean; default: TRUE.
#'              Whether rows that contains only null value should be removed
#'
#' @return  a data.frame which contains only the columns that contains the tag
#' @export
#'
extractWithinSampleAbundance<-extractWithTag(tag="WithInSampleAbundance.")

#' Homogenise the names of the Uhlen et al. data with the other datasets
#' (in particular with GTEx)
#'
#' @param DF a Uhlen expression (numeric) data.frame
#'
#' @return a data.frame with homogenised columns names
#' @export
#'
fixNamesUhlen<-function(DF){
  colnames(DF)<-gsub('Endometrium','Uterus',colnames(DF))
  colnames(DF)<-gsub('Cerebral.cortex','Cortex',colnames(DF))
  colnames(DF)<-gsub('Esophagus','Oesophagus',colnames(DF))
  colnames(DF)<-gsub('[.]',' ',colnames(DF))
  colnames(DF)<-gsub('Gallbladder','Gall bladder',colnames(DF))
  return(DF)
}

#' Homogenise the names of the GTEx data with the other datasets
#' (in particular with Uhlen et al. )
#'
#'
#' @param DF a GTEx expression (numeric) data.frame
#'
#' @return a data.frame with homogenised columns names
#' @export
#'
fixNamesGtex<-function(DF){
  colnames(DF)<-gsub('Esophagus','Oesophagus',colnames(DF))
  colnames(DF)<-gsub('Frontal.cortex','Frontalcortex',colnames(DF))
  colnames(DF)<-gsub('[.]',' ',colnames(DF))
  colnames(DF)<-gsub('Gallbladder','Gall bladder',colnames(DF))
  return(DF)
}
