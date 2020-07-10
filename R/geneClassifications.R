#' Give the list of rows for which the values are only once greater than the threshold
#'
#' @param DF a numeric data.frame
#' @param threshold numeric; default 0
#' @param type "gt" for greater than or anything else for greater than or equal to
#'
#' @return list of row names
#' @export
#'
protSpe<-function(DF,threshold=0,type='gt'){

  end<-ncol(DF)
  DF["name"]<-rownames(DF)

  if(type=='gt'){
    list.unique<-apply(DF,1, function(x) {
      y<-x[1:end]
      if ( sum( as.numeric(x[1:end]) > threshold) == 1) {
        return(setNames(names(y)[which(as.numeric(y)>0)],x['name']))
      }else{
        return(setNames(NA,x[c(length(x))]))
      }
    })
    return(list.unique)
  }else{
    list.unique<-apply(DF,1, function(x) {
      y<-x[1:end]
      if (sum(as.numeric(x[1:end]) %>=% threshold)==1) {
        return(setNames(names(y)[which(as.numeric(y)>0)],x['name']))
      }else{
        return(setNames(NA,x[c(length(x))]))
      }
    })
    return(list.unique)
  }
}

#' Classify the genes in function of their expression to be used with clustermq. Inspired from Uhlén et al. classification.
#' @param DF numeric data.frame
#' @param detect numeric; default:0. Threshold to consider the gene expressed
#' @param groupNB parameter
#' @param HighThreshold numeric. Threshold above which the expression is considered high.
#'                      default: 10 as it is in Uhlén et al. papers
#' @param foldLow numeric. Expression fold to be considered as 'TissueSpecific'. Default:5
#' @param foldHigh numeric. Expression fold to be considered as 'TissueSpecific.High'. Default: 50/
#'
#' @return vector of genes names
#' @export
geneClassify.Q<-function(DF,detect=0,groupNB,HighThreshold=10,
                         foldLow=5,foldHigh=50){
  if(!base::missing(groupNB)) groupNB=round(nrow(DF)/4)
  bSup<-nrow(DF)
  GeneID<-colnames(DF)#geneID is the only column

  DF<-data.frame(DF)

  if(detect==0){
    nbSums<-colSums(DF[1,]>detect)
  }else{
    nbSums<-colSums(DF[1,]>=detect)
  }
  if(nbSums==0) return(setNames("undetected",GeneID))
  if(nbSums==1) return(setNames("TissueSpe",GeneID))

  means<-colMeans(DF[1,])

  if(max(DF[1,])/means>=foldHigh) return(setNames('TissueSpecific.High',GeneID))
  if(max(DF[1,])/means>=foldLow) return (setNames('TissueSpecific',GeneID))
  if(any(DF[1,]>5*means)){
    if(max(sort(DF[1,])[1:(bSup-groupNB)])<5*means)
      return(setNames('GroupEnriched',GeneID))
  }
  if(nbSums==bSup) {
    if(min(DF[1,])>=HighThreshold){
      return(setNames('Ubiquitous.High',GeneID))
    }else{
      return(setNames('Ubiquitous.Low',GeneID))
    }
  }else{
    if(max(sort(DF[1,],decreasing=TRUE)[1:nbSums])<HighThreshold){
      return(setNames('Mixed.Low',GeneID))
    }else{
      return(setNames('Mixed.High',GeneID))
    }
  }

}

#' Classify the genes in function of their expression. Inspired from Uhlén et al. classification.
#'
#' @param DF numeric data.frame
#' @param detect numeric; default:0. Threshold to consider the gene expressed
#' @param groupNB parameter
#' @param HighThreshold numeric. Threshold above which the expression is considered high.
#'                      default: 10 as it is in Uhlén et al. papers
#' @param foldLow numeric. Expression fold to be considered as 'TissueSpecific'. Default:5
#' @param foldHigh numeric. Expression fold to be considered as 'TissueSpecific.High'. Default: 50/
#'
#' @return vector of genes names
#' @export
#'
geneClassify<-function(DF,detect=0,groupNB,HighThreshold=10,
                       foldLow=5,foldHigh=50){
  if(!base::missing(groupNB)) groupNB=round(ncol(DF)/4)
  bSup<-ncol(DF)
  DFtemp<-parallel::mclapply(rownames(DF),function(x){
    sDF<-DF[x,]
    if(detect==0){
      nbSums<-rowSums(sDF>detect)
    }else{
      nbSums<-rowSums(sDF%>=% detect)
    }

    if(nbSums==0) return('Undetected')
    if(nbSums==1) {
      if(max(sDF)>foldHigh) {
        return('TissueSpecific.High')
      }else{
        return('TissueSpecific')
      }
    }
    means<-rowMeans(sDF)
    if(max(sDF)/means>=foldHigh) return('TissueSpecific.High')
    if(max(sDF)/means>=foldLow) return ('TissueSpecific')
    if(any(sDF>5*means)){
      if(max(sort(sDF)[1:(bSup-groupNB)])<5*means)
        return('GroupEnriched')
    }
    if(nbSums==bSup) {
      if(min(sDF)>=HighThreshold){
        return('Ubiquitous.High')
      }else{
        return('Ubiquitous.Low')
      }
    }else{
      if(max(sort(sDF,decreasing=TRUE)[1:nbSums])<HighThreshold){
        return('Mixed.Low')
      }else{
        return('Mixed.High')
      }
    }

  })
  return(c(unlist(setNames(DFtemp,rownames(DF)))))
}

#' Rank for each tissue (column) the genes (rows)
#' based on their relative specificity to that tissue
#'
#' @param DF numeric data.frame
#' @param detect numeric. Threshold for detection; default: 0
#'
#' @return numeric data.frame
#' @export
#'
rank.Atlas<-function(DF,detect=0){
  nameRow<-rownames(DF)
  size<-dim(DF)[2]-1
  resDF<-lapply(colnames(DF),function(x){
    #vector DF[,x] is going to be compared with the corresponding row of DF[,names(DF)!=x]
    listrow<-lapply(nameRow,function(y) {
      if (detect){#a threshold has been given,
        #only the genes higher than this threshold are going to be used
        if(DF[y,x]>=detect){
          #	if (DF[y,x]==max(DF[y,]) #we are only interrested
          value<-DF[y,x]/max(DF[y,names(DF)!=x]) #min fold that the current value has compared to the next max

        }else{#not detected
          value<--1
        }
      }else{
        value<-DF[y,x]/max(DF[y,names(DF)!=x]) #min fold that the current value has compared to the next max
      }
      setNames(value,y)

    }
    )
    resCond<-unlist(c(listrow))
    dd.out<-as.data.frame(resCond)
    colnames(dd.out)<-x
    dd.out
  }
  )
  res<-Reduce(cbind,resDF)
  return(res)
}


#' Rank expression of the genes within each tissue.
#'
#' @param DF numeric data.frame
#'
#' @return numeric data.frame
#' @export
#'
ebiRank<-function(DF){
  n<-ncol(DF)
  finalRes<-sapply(rownames(DF),function(x){
    y<-DF[x,]
    tmp<-sort(y)
    MAX<-tmp[n]
    MAX2<-tmp[n-1]
    y[]<-sapply(y,function(x) ifelse(x<MAX,x/MAX,x/MAX2))
    return(y)
  })
  return(t(finalRes))
}


#' Rank gene expression based on their specificity to each tissue
#'
#' @param DF numeric data.frame
#' @param ties.method how to handle ties. eg: 'min' or 'max'
#' @param empty how to handle missing data or bad score
#' @param order logical; default: TRUE. Whether the result should be order
#' @param log2 logical; default:FALSE. Whether the data.frame should first be transformed
#'
#' @return numeric data.frame
#' @export
rankSpe<-function(DF,ties.method='min',empty='rank',order=TRUE,log2=FALSE){
  ##the function gives back a data.frame or matrix with the same dimensions
  ##and same column and row names

  DF.raw<-DF
  ord.orig.row<-rownames(DF)

  ##First, prot/genes detected in one condition only or not at all.
  uniqCond<-rownames(DF[rowSums(DF)<2,])

  if(log2) DF<-log2(DF+1)
  DFres<-data.frame(lapply(colnames(DF),function(Colx){
    #for each of the condition:
    #initialise tempRes
    tempRes<-NULL
    #the ones that have a specific expression
    uniq<-setNames(DF.raw[uniqCond,Colx],uniqCond)
    #the following one will be scored by their expression data
    uniqExp<-uniq[uniq>0]
    if(exists("uniqExp")){
      uniqExp<-base::rank(-uniqExp,ties.method=ties.method)
      tempRes<-c(tempRes,uniqExp)
    }
    #these ones will have NA or a very bad score at the end.
    nillExp<-uniq[!names(uniq) %in% names(uniqExp)]
    if(exists('nillExp')){
      if(empty=='rank'){
        nillExp[nillExp==0]<-nrow(DF)-length(nillExp)+1
      }else{
        nillExp[nillExp==0]<-NA
      }
      tempRes<-c(tempRes,nillExp)
    }
    #the others
    others<-setNames(DF[,Colx],rownames(DF))
    others<-others[!names(others) %in% uniqCond]
    others<-others+1
    others<-others/rowSums(DF[names(others),]+1)
    others<-base::rank(-others,ties.method=ties.method)
    if(exists('uniqExp')) others<-others+length(uniqExp)

    tempRes<-c(tempRes,others)
    #final check
    if(length(tempRes)<length(ord.orig.row)){
      missingVal<-ord.orig.row[!ord.orig.row %in% names(tempRes)]
      missingVal<-setNames(rep(NA,length(missingVal)),missingVal)
      tempRes<-c(tempRes,missingVal)
    }
    return(tempRes[ord.orig.row])
  })#for each of the condition, each protein/gene has been ranked and then put back in the original order
  )#dataframe

  colnames(DFres)<-colnames(DF)
  nb<-ncol(DFres)
  DFres$score.mean<-rowMeans(DFres[,1:nb])
  DFres$score.sd<-rowSD(DFres[,1:nb])
  if(order) DFres<-DFres[order(DFres$score.sd),]

  return(DFres)
}


#' Give the gene specificity in each tissue as a ratio to its mean across all tissue
#' @description to avoid division by 0, 1 is added to all values.
#'
#' @param DF numeric data.frame
#'
#' @return numeric data.frame
#' @export
#'
ratioSpe.1<-function(DF){
  return((DF+1)/rowMeans(DF+1))
}

# a priori used in the thesis for the transcriptomics chapter
#' Compute the tissue specificity of the genes to each tissue as a ratio
#' to the sum of the gene expression across all tissues.
#'
#' @param DF numeric data.frame
#'
#' @return numeric data.frame
#' @export
#'
ratioSpe<-function(DF){
  finalRes<-Reduce(rbind,lapply(rownames(DF),function(x){
    SUM<-rowMeans(DF[x,])
    if(SUM>0){
      res<-DF[x,]/SUM
    }else{
      res<-rep(-1,ncol(DF))
    }}))
  finalRes<-as.data.frame(finalRes)
  return(finalRes)
}



# Classification inspired by Faberger et al. and Uhlén et al. --------

#' List all the genes (rows) for which no expression was detected at a given threshold
#'
#' @param DF numeric data.frame
#' @param detect numeric; threshold of detection; default: 1
#'
#' @return vector of characters (gene names)
#' @export
#'
undetected<-function(DF,detect=1){
  rownames(DF[apply(DF, 1, function(x) {all(x <detect)}),])
}


#' Give the names of all the genes (rows) that are strictly expressed in every tissue (column) above a given threshold
#'
#' @param DF numeric data.frame
#' @param threshold numeric; default: 0. Minimal threshold above which the gene has to be observed
#'                  to be considered as expressed.
#'
#' @return List of the names of the rows complying with the threshold in all tissues (columns)
#' @export
#'
ubi.strict<-function(DF,threshold=0){
  end<-dim(DF)[2]
  DF["name"]<-rownames(DF)
  list.high<-apply(DF,1, function(x) {
    if (all(as.numeric(x[1:end])>threshold)) {
      x[c(length(x))] #send back the rownames that is stored in the last column
    }
  })
  return(unname(Reduce(c,list.high)))
}

#' Give the genes (rows) that are ubiquitously highly expressed across all tissues (columns)
#'
#' @param DF numeric data.frame
#' @param threshold minimal threshold; default: 10
#'
#' @return vector with the names of the rows that check the condition
#' @export
#'
ubi.high<-function(DF,threshold=10.0){
  end<-dim(DF)[2]
  DF["name"]<-rownames(DF)
  list.high<-apply(DF,1, function(x) {
    if (all(as.numeric(x[1:end])>=threshold)) {
      x[c(length(x))]
    }
  })
  if(length(list.high)>0) {
    return(unname(Reduce(c,list.high)))
  }else{
    return(character(0))
  }
}


#' Give the genes (rows) that are ubiquitously lowly expressed across all tissues (columns)
#'
#' @param DF a numeric data.frame
#' @param threshold maximal threshold of expression; default: 10
#' @param detect minimal level of expression to be considered as expressed; default:1
#'
#' @return vector with the names of the rows that are
#' @export
#'
ubi.low<-function(DF,threshold=10,detect=1){
  end<-base::ncol(DF)
  DF["name"]<-rownames(DF)
  list.low<-apply(DF,1, function(x) {
    if (all(as.numeric(x[1:end]) <threshold &as.numeric(x[1:end]) >detect)) {
      x[c(length(x))]
    }
  })
  if(length(list.low)>0){
    return(unname(Reduce(c,list.low)))
  }else{
    return(character(0))
  }
}

#' Provide the list of mixed high genes
#'
#' @description For each gene (rows) the number of samples/tissues (columns)
#' in which they are expressed is given as the content of the named vector.
#' Genes have to be consistantly highly expressed (ie expressed above a given threshold)
#' when they are detected (ie genes have to be expressed above another given threshold).
#'
#' @param DF numeric data.frame
#' @param threshold numeric. Minimal level of expression when detected; default: 10.0
#' @param detect numeric. Minimal level to be considered as detected; default: 1
#'
#' @return a named vector with the number of samples the genes has been detected as mixed high
#' @export
#'
mixed.high<-function(DF,threshold=10.0,detect=1){
  end<-base::ncol(DF)
  DF["name"]<-rownames(DF)
  list.mixhigh<-apply(DF,1, function(x) {
    if (all(as.numeric(x[1:end])<detect|as.numeric(x[1:end])>=threshold)) {
      setNames(sum(as.numeric(x[1:end])>=threshold),x[c(length(x))])
    }
  })
  result<-Reduce(c,list.mixhigh)
  return(result[result!=end])
}


#' Provide the list of genes that have a high expression when they are expressed
#'
#' @param DF numeric data.frame
#' @param threshold numeric. Minimal level of expression when detected; default: 10.0
#' @param detect numeric. Minimal level to be considered as detected; default: 1
#' @param clean logical; default: TRUE. Wheter the genes that aren't mixed high
#'              should be removed from the final result.
#'
#' @return named vector. Names of the genes that have a mixed high expression.
#' @export
#'
mixed.high2<-function(DF,threshold=10.0,detect=1,clean=TRUE){
  finalRes<-sapply(rownames(DF),function(x){
    if(all(DF[x,]<detect|DF[x,]%>=%threshold)){
      return(paste(colnames(DF)[which(DF[x,]>threshold)],collapse = ';'))
    }else{
      return('')
    }
  })
  if(clean) finalRes<-finalRes[finalRes!='']
}


#' Provide the list of mixed low genes
#'
#' @description For each gene (rows) the number of samples/tissues (columns)
#' in which they are expressed is given as the content of the named vector.
#' Genes have to be consistantly lowly expressed (ie expressed lower a given threshold)
#' when they are detected (ie genes have to be expressed above another given threshold).
#'
#' @param DF a numeric data.frame
#' @param threshold Maximal level of expression when detected; default: 10.0
#' @param detect numeric. Minimal level to be considered as detected; default: 1
#'
#' @return a named vector with the number of samples the genes has been detected as mixed low
#' @export
#'
mixed.low<-function(DF,threshold=10.0,detect=1){
  end<-base::ncol(DF)
  DF["name"]<-rownames(DF)
  list.mixlow<-apply(DF,1, function(x) {
    if (all(as.numeric(x[1:end])<detect|as.numeric(x[1:end])<threshold)) {
      setNames(sum(as.numeric(x[1:end])>detect),x[c(length(x))])
      #return a vector with the number of samples the genes has been detected
      #the gene has to been expressed above the threshold whenever detected
    }
  })
  result<-Reduce(c,list.mixlow)
  return(result[result!=end])
}


#' Provide the list of genes that have a low expression when they are expressed
#'
#' @param DF a numeric data.frame
#' @param threshold Maximal level of expression when detected; default: 10.0
#' @param detect numeric. Minimal level to be considered as detected; default: 1
#' @param clean logical; default: TRUE. Wheter the genes that aren't mixed low
#'              should be removed from the final result.
#'
#' @return named vector. Names of the genes that have a mixed low expression.
#' @export
#'
mixed.low2<-function(DF,threshold=10.0,detect=1,clean=TRUE){
  finalRes<-sapply(rownames(DF),function(x){
    if(all(DF[x,]<detect|DF[x,]<threshold)){
      return(paste(colnames(DF)[which(DF[x,]<threshold)],collapse = ';'))
    }else{
      return('')
    }
  })
  if(clean) finalRes<-finalRes[finalRes!='']
}


#' Annotate the genes that the expression is enriched of enhanced for a specific tissue
#'
#' @param DF numeric data.frame
#' @param detect numeric; default:0. Minimal level of expression to be considered detected
#' @param factors logical. Whether the outputed classification should be characters or factors; default: FALSE
#'
#' @return data.frame
#' @export
#'
relative.GXenrichement<-function(DF,detect=0,factors=FALSE){
  nameRow<-rownames(DF)
  size<-dim(DF)[2]-1
  if (detect) { # a threshold has been given
    resDF<-parallel::mclapply(colnames(DF),function(x){
      #vector DF[,x] is going to be compared with the corresponding row of DF[,names(DF)!=x]
      listrow<-parallel::mclapply(nameRow,function(y) {
        #only the genes higher than the threshold "detect" are going to be used
        MEAN<-rowMeans(DF[y,])
        if(DF[y,x]>=detect){
          maxR<-max(DF[y,names(DF)!=x])

          if (DF[y,x]>5*maxR){
            if(DF[y,x]>50*maxR) value<-">50x" else value<-">5x"
          }else{
            if(DF[y,x]>5*MEAN){
              value<-'Tissue Enhanced'
            }else{
              if(!DF[y,x]*5<maxR){
                nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
                value<-paste(nbTrue,size,sep='/')
              }else{
                value<-0
              }
            }
          }
        }else{#not detected
          value<--1
        }
        setNames(value,y)
      }
      )
      resCond<-unlist(c(listrow))
      dd.out<-as.data.frame(resCond)
      colnames(dd.out)<-x
      dd.out
    }
    )
    res<-Reduce(cbind,resDF)
  }else{ #no threshold given at the function call
    resDF<-parallel::mclapply(colnames(DF),function(x){
      #vector DF[,x] is going to be compared with the corresponding row of DF[,names(DF)!=x]
      listrow<-parallel::mclapply(nameRow,function(y) {
        #only the genes higher than the threshold "detect" are going to be used
        maxR<-max(DF[y,names(DF)!=x])
        MEAN<-rowMeans(DF[y,])

        if (DF[y,x]>5*maxR){
          if(DF[y,x]>50*maxR) value<-">50x" else value<-">5x"
        }else{
          if(DF[y,x]>5*MEAN){
            value<-'Tissue Enhanced'
          }else{
            if(!DF[y,x]*5<maxR){
              nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
              value<-paste(nbTrue,size,sep='/')
            }else{
              value<-0
            }
          }
        }
        setNames(value,y)
      }
      )
      resCond<-unlist(c(listrow))
      dd.out<-as.data.frame(resCond)
      colnames(dd.out)<-x
      dd.out
    }
    )

    res<-Reduce(cbind,resDF)
  }

  if(!factors){
    return(factorDF2string(res))
  }else{
    return(res)
  }
}

#' Another function that annotate the genes that the expression is enriched of enhanced for a specific tissue
#'
#' @param DF numeric data.frame
#' @param detect numeric; default:0. Minimal level of expression to be considered detected
#' @param factors logical. Whether the outputed classification should be characters or factors; default: FALSE
#'
#' @return data.frame
#' @export
#'
relative.GXenrichement2<-function(DF,detect=0,factors=FALSE){
  nameRow<-rownames(DF)
  size<-dim(DF)[2]-1
  resDF<-parallel::mclapply(colnames(DF),function(x){
    #vector DF[,x] is going to be compared with the corresponding row of DF[,names(DF)!=x]
    listrow<-lapply(nameRow,function(y) {
      if (detect){#a threshold has been given,
        #only the genes higher than this threshold are going to be used
        if(DF[y,x]>=detect){
          nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
          if (nbTrue==size){
            if (all(DF[y,names(DF)!=x]*50<DF[y,x])) value<-">50x" else value<-">5x"
          }else{
            if(nbTrue>1) value<-paste(nbTrue,size,sep='/') else value<-0
          }
        }else{#not detected
          value<--1
        }
      }else{
        nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
        if (nbTrue==size){
          if (all(DF[y,names(DF)!=x]*50<DF[y,x])) value<-">50x" else value<-">5x"
        }else{
          if(nbTrue>1) value<-paste(nbTrue,size,sep='/') else value<-0
        }
      }
      setNames(value,y)
    }
    )
    resCond<-unlist(c(listrow))
    dd.out<-as.data.frame(resCond)
    colnames(dd.out)<-x
    dd.out
  }
  )
  res<-Reduce(cbind,resDF)
  if(!factors){
    return(factorDF2string(res))
  }else{
    return(res)
  }
}



#' Another function that annotates the genes that the expression is enriched of enhanced for a specific tissue,
#' but the output is numeric
#'
#' @param DF numeric data.frame
#' @param detect numeric; default:0. Minimal level of expression to be considered detected
#'
#'
#' @return numeric data.frame
#' @export
#'
relnum.GXenrichement<-function(DF,detect=0){
  nameRow<-rownames(DF)
  size<-dim(DF)[2]-1
  resDF<-lapply(colnames(DF),function(x){
    #vector DF[,x] is going to be compared with the corresponding row of DF[,names(DF)!=x]
    listrow<-lapply(nameRow,function(y) {
      if (detect){#a threshold has been given,
        #only the genes higher than this threshold are going to be used
        if(DF[y,x]>=detect){
          nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
          if (nbTrue==size){
            if (all(DF[y,names(DF)!=x]*50<DF[y,x])) value<-50*size else value<-5*size
          }else{
            if(nbTrue>1) value<-nbTrue else value<-0
          }
        }else{#not detected
          value<--1
        }
      }else{
        nbTrue<-sum(DF[y,names(DF)!=x]*5<DF[y,x])
        if (nbTrue==size){
          if (all(DF[y,names(DF)!=x]*50<DF[y,x])) value<-50*size else value<-5*size
        }else{
          if(nbTrue>1) value<-nbTrue else value<-0
        }
      }
      setNames(value,y)
    }
    )
    resCond<-unlist(c(listrow))
    dd.out<-as.data.frame(resCond)
    colnames(dd.out)<-x
    dd.out
  }
  )
  res<-Reduce(cbind,resDF)
  return(res)
}



#' Give the list of the genes whose expression is enriched for a tissue at the given threshold.
#'
#' @param DF numeric data.frame
#' @param fold numeric. default: 5. Above which fold change a gene is considered as enriched for a tissue
#' @param detect numeric; default:0. Minimal level of expression to be considered detected
#' @param threshold numeric; default:1. Minimal level of expression to be considered as possibly translated as protein
#' @param verbose logical; default:TRUE. Print to screen the number of genes that are enriched and aren't
#' @param clean logical; default: TRUE. Remove the genes that haven't show any enrichment
#'
#' @return named vector. Names give the gene ID (from the rownames)
#'         and the content is the tissue where the enrichment has been observed
#' @export
#'
enrichedTissue<-function(DF,fold=5,detect=FALSE,threshold=1,verbose=TRUE,clean=TRUE){
  if(detect){
    if(threshold>0)
      DF<-strip(DF,'ge',threshold)
    else
      DF<-strip(DF,'gt',0)
  }
  n<-ncol(DF)
  finalRes<-sapply(rownames(DF),function(x){
    tmp<-sort(DF[x,])
    if(tmp[n]>fold*tmp[n-1])
      res<-names(tmp[n])
    else
      res<-''
  })
  if(verbose) print(summary(as.factor(finalRes)))
  if(clean) finalRes<-finalRes[finalRes!='']

  return(finalRes)
}


#' Give the list of the genes whose expression is enhanced for a tissue at the given threshold.
#'
#' @param DF numeric data.frame
#' @param fold numeric; default: 5. Above which fold change a gene is considered as enhanced for a tissue
#' @param detect numeric; default:0. Minimal level of expression to be considered detected
#' @param threshold numeric; default:1. Minimal level of expression to be considered as possibly translated as protein
#' @param clean logical; default: TRUE. Remove the genes that haven't show any enhancement
#'
#' @return named vector. Names give the gene ID (from the rownames)
#'         and the content is the tissue where the enhancement has been observed
#' @export
#'
enhancedTissue<-function(DF,fold=5,detect=FALSE,threshold=1,clean=TRUE){
  if(detect){
    if(threshold>0)
      DF<-strip(DF,'ge',threshold)
    else
      DF<-strip(DF,'gt',0)
  }

  means<-rowMeans(DF)

  finalRes<-sapply(rownames(DF),function(x){
    return(paste(colnames(DF)[which(DF[x,]> fold*means[x])],collapse=';'))
  })

  if(clean) finalRes<-finalRes[finalRes!='']

  return(finalRes)

}

#' Create one data.frame with the most relevant information for each gene
#' from previous gene classification ("Uhlén-like) analysis
#'
#' @param notExpressed object, or filename
#' @param threshold level of expression
#' @param NotDetected object, or filename
#' @param Enriched object, or filename
#' @param Enhanced object, or filename
#' @param ubiHigh object, or filename
#' @param ubiLow object, or filename
#' @param mixedHigh object, or filename for not expressed genes
#' @param mixedLow object, or filename for not expressed genes
#' @param outfilename filename for the resulting data.frame to be saved with
#' @param gene.mapID Map of the gene names with their Ensembl ids
#' @param precomputed logical, whether to search for recorded data
#'
#' @return a data.frame with all the genes annotated with their "Uhlén-like" classification.
#'         result can be saved directly to a file simultaneously
#' @export
#'
sortUhlenClassInOne<-function(notExpressed,threshold,NotDetected,Enriched,
                              Enhanced,ubiHigh,ubiLow,mixedHigh,mixedLow,
                              outfilename,gene.mapID,precomputed=FALSE){

  addEmptyDF<-function(DF,gene,Class){
    return(rbind(DF,
                data.frame(GeneID=gene,Name=NA,Class=Class,
                           Tissues=NA,stringsAsFactors = FALSE)))
  }

  DF<-data.frame(GeneID=character(0),Name=character(0),
                 Class=character(0),Tissues=character(0),
                 stringsAsFactors = FALSE)

  if(precomputed){
    if(isTRUE(grep('.RData',notExpressed[1]))) notExpressed<-Load(notExpressed)
    if(isTRUE(grep('.RData',notDetected[1]))) NotDetected<-Load(NotDetected)
    if(isTRUE(grep('.RData',Enriched[1]))) Enriched<-Load(Enriched)
    if(isTRUE(grep('.RData',Enhanced[1]))) Enhanced<-Load(Enhanced)
    if(isTRUE(grep('.RData',ubiHigh[1]))) ubiHigh<-Load(ubiHigh)
    if(isTRUE(grep('.RData',ubiLow[1]))) ubiLow<-Load(ubiLow)
    if(isTRUE(grep('.RData',mixedHigh[1]))) mixedHigh<-Load(mixedHigh)
    if(isTRUE(grep('.RData',mixedLow[1]))) mixedLow<-Load(mixedLow)
  }

  if(!missing(Enriched)){
    if(length(Enriched)){
      DF<-data.frame(GeneID=names(Enriched),
                     Name=gene.mapID[names(Enriched)],
                     Class='Tissue Enriched',
                     Tissues=Enriched,
                     stringsAsFactors = FALSE)
    }else{
      DF<-addEmptyDF(DF,gene="NoEnriched",Class="Tissue Enriched")
    }
  }

  if(!missing(Enhanced)){
    if(length(Enhanced)){
      Enhanced<-Enhanced[!(names(Enhanced) %in% DF$GeneID)]
      if(length(Enhanced)) DF<-rbind(DF,
                                     data.frame(GeneID=names(Enhanced),
                                     Name=gene.mapID[names(Enhanced)],
                                     Class='T/G Enhanced',
                                     Tissues=Enhanced,
                                     stringsAsFactors = FALSE))
    }
    if(!length(Enhanced)){
      DF<-addEmptyDF(DF,gene='NoEnhanced',Class='T/G Enhanced')
    }
  }

   if(!missing(ubiHigh)){
     if(length(ubiHigh)){
       ubiHigh<-ubiHigh[!(ubiHigh %in%  DF$GeneID)]
       if(length(ubiHigh)) DF<-rbind(DF,
                                     data.frame(GeneID=ubiHigh,
                                     Name=gene.mapID[ubiHigh],
                                     Class='Highly expressed ubiquitously',
                                     Tissues='All',
                                     stringsAsFactors = FALSE))

     }
     if(!length(ubiHigh)){
       DF<-addEmptyDF(DF,gene='NoUbiHigh',Class='Highly expressed ubiquitously')
     }
   }

  if(!missing(ubiLow)){
    if(length(ubiLow)){
      ubiLow<-ubiLow[!(ubiLow %in% DF$GeneID)]
      if(length(ubiLow)) DF<-rbind(DF,
                                   data.frame(GeneID=ubiLow,
                                   Name=gene.mapID[ubiLow],
                                   Class='Lowly expressed ubiquitously',
                                   Tissues='All',
                                   stringsAsFactors = FALSE))

    }
    if(!length(ubiLow)){
      DF<-addEmptyDF(DF,gene='NoUbiLow', Class='Lowly expressed ubiquitously')
    }
  }

  if(!missing(mixedHigh)){
    if(length(mixedHigh)){
      mixedHigh<-mixedHigh[!(names(mixedHigh) %in% DF$GeneID)]
      if(length(mixedHigh)) DF<-rbind(DF,
                                       data.frame(GeneID=names(mixedHigh),
                                       Name=gene.mapID[names(mixedHigh)],
                                       Class='Highly expressed when detected',
                                       Tissues=mixedHigh,
                                       stringsAsFactors = FALSE))
    }
    if(!length(mixedHigh)){
      DF<-addEmptyDF(DF,gene='NoMixedHigh',Class='Highly expressed when detected')
    }
  }

  if(!missing(mixedLow)){
    if(length(mixedLow)){
      mixedLow<-mixedLow[!(names(mixedLow) %in% DF$GeneID)]
      if(length(mixedLow)) DF<-rbind(DF,
                                     data.frame(GeneID=names(mixedLow),
                                     Name=gene.mapID[names(mixedLow)],
                                     Class='Lowly expressed ubiquitously',
                                     Tissues=mixedLow,
                                     stringsAsFactors = FALSE))

    }
    if(!length(mixedLow)){
      DF<-addEmptyDF(DF,gene='NoMixedLow',Class="Lowly expressed ubiquitously")
    }
  }

  if(threshold>0){#otherwise the same as not detected
    if(!missing(notExpressed)){
      if(length(notExpressed)){
        notExpressed<-notExpressed[!(notExpressed %in%  DF$GeneID)]
        if(length(notExpressed)) DF<-rbind(DF,
                                           data.frame(GeneID=notExpressed,
                                           Name=gene.mapID[notExpressed],
                                           Class=paste('Not expressed at',threshold),
                                           Tissues='None',stringsAsFactors = FALSE))

      }
      if(!length(notExpressed)){
        DF<-addEmptyDF(DF, gene='NogeneNotExpressed',Class=paste('Not expressed at',threshold))
      }
    }
  }

  if(!missing(NotDetected)){
    if(length(NotDetected)){
      NotDetected<-NotDetected[!(NotDetected %in% DF$GeneID)]
      if(length(NotDetected)) DF<-rbind(DF,
                                        data.frame(GeneID=NotDetected,
                                                   Name=gene.mapID[NotDetected],
                                                   Class='Not detected',
                                                   Tissues='None',stringsAsFactors = FALSE))
    }
    if(!length(notDetected)){
      DF<-rbind(DF,
                data.frame(GeneID='NoneNotDetected',
                           Name=NA,Class='Not detected',
                           Tissues=NA,stringsAsFactors = FALSE))
    }
  }

  if(any(!is.na(DF$GeneID))){
    rownames(DF)<-DF$GeneID
  }

  if(!missing(outfilename)) saveToFile(DF,filename = outfilename)
  return(DF)

}


#' Give the list of expressed genes in a study in function of the selected type of comparison
#'
#' @param DF numeric data.frame
#' @param threshold numeric, default: 0. Threshold above which a gene is considered as expressed
#' @param type type of comparison to be considered.
#'             'gt': expression should be strictly greater than threshold.
#'             'ge': expression should be greater than or equal to threshold.
#'             'eq': expression should be equal to threshold.
#'             'lt': expression should strictly be lesser than threshold
#'             'le': expression should be lesser than or equal to threshold
#'             'notDetected': No expression is detected above or at threshold
#'             'ubi': expression is detected everywhere
#'             'unique': expression is detected solely in one condition
#'             'mixedH': expression is detected above threshold in 'nb' tissues
#'             'mixedl': expression is detected below threshold in 'nb' tissues
#' @param out result format;
#'            'rownames': names of the relevant genes
#'            'data.frame': data.frame with the relevant genes
#'            'df': alias of 'data.frame'
#'            'index': indices of the rows of the relevant genes
#'            'ind': alias of 'index'
#' @param nb Number of tissues that one should be expressed to be considered for mixed expression
#'
#' @return several types of objects are available. Selection through "out"
#' @export
#'
ListOfExpressedGenes<-function(DF,threshold=0,out='rownames',type='gt',nb){
  #Take a data.frame and send back rownames or rows that have an expression above a certain threshold
  #The DF has to be only numeric or integer
  #try(DF <- as.numeric(DF), stop("first element need to be all numeric"))

  DF1<-'No result for this threshold' #to deal with the output
  # Are only selected the rows which follow the specified rule
  switch(type,
         'gt'={ DF1<-DF[rowSums(DF > threshold) !=0,] },
         'ge'={ DF1<-DF[rowSums(DF %>=% threshold) !=0,] },
         'eq'={ DF1<-DF[rowSums(DF == threshold) !=0,] },
         'lt'={ DF1<-DF[rowSums(DF < threshold) !=0,] },
         'le'={ DF1<-DF[rowSums(DF %<=% threshold) !=0,] },
         'notDetected'={ DF1<-DF[rowSums(DF %<=% threshold) == ncol(DF) ,] },
         'ubi'={ DF1<-DF[rowSums(DF %>=% threshold) == ncol(DF) ,] },
         'unique'={ DF1<-DF[rowSums(DF %>=% threshold) == 1, ] }, #specific case of `mixed`
         'mixedH'={ DF1<-DF[rowSums(DF %>=% threshold) == nb, ] },
         'mixedl'={ DF1<-DF[rowSums(DF %<=% threshold) == nb, ] }
  )
  if(class(DF1)[1]=='character') out<-'noResult' #bad! should be changed

  switch(out,
         'rownames'=return(rownames(DF1)),
         'data.frame'=return(DF[rownames(DF1),]),
         'df'=return(DF[rownames(DF1),]),
         'index'=return(rownames(DF) %in% rownames(DF1)),
         'ind'=return(rownames(DF) %in% rownames(DF1)),
         'noResult'=return(DF1)
  )
}

