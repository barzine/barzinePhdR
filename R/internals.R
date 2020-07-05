prep<-function(DF,threshold,label,midLabel,endLabel,out='DF'){

  DF<-extractSpe(DF,threshold,verbose=TRUE,format='lg')

  if(!missing(midLabel)&&!missing(endLabel)){
    if(threshold==0){
      DF$Label<-paste0(label,' ',midLabel,'> ',threshold,' ',endLabel)
    }else{
      DF$Label<-paste0(label,' ',midLabel,"\u2265 ",threshold,' ',endLabel)
    }
  }else{
    DF$Label<-label
  }
  switch(out,
         'DF'   = return(DF),
         'list' = return(list(DF,unique(DF$Label)))
         )
}

prep2<-function(d.uvar,DF,decreas){
  sortD<-setNames(rep(0,length(d.uvar)),d.uvar)
  sortD<-sapply(names(sortD),function(x){
    nrow(DF[DF$variable==x,])})
  sortD<-sort(sortD,decreasing=decreas)
  return(sortD)
}



sharedBreadth_firstLast<-function(...,a,b,typeR='vec',colNameDF='nb.tissues',indexCol,threshold){
  if(missing('threshold')) threshold=0
  if(missing('indexCol')) indexCol<-c(1,2,min(ncol(a),ncol(b))-1)
  unbalanced<-FALSE
  if(ncol(a)!=ncol(b)) {
    unbalanced=TRUE
    newIndexCol<-setdiff(indexCol,min(ncol(a),ncol(b)))
  }

  a$shared<-ifelse(rownames(a) %in% rownames(b),'Shared','Unshared')
  a$shared<-unlist(lapply(rownames(a),function(x){
    if(!unbalanced){
      if(a[x,'shared']=='Shared'){
        if(a[x,colNameDF] %in% indexCol & a[x,colNameDF]==b[x,colNameDF]){
          return('Identical')
        }else{
          return('Different')
        }
      }else{
        return('Unshared')
      }
    }else{
      if(a[x,'shared']=='shared'){
        if(a[x,colNameDF] %in% newIndexCol & a[x,colNameDF]==b[x,colNameDF]){
          return('Identical')
        }else{
          return('Different')
        }
        if(a[x,colNameDF]==indexCol[length(indexCol)]){
          if(b[x,colNameDF]==indexCol[length(indexCol)]){
            return('Identical')
          }else{
            if(b[x,colNameDF]<indexCol[length(indexCol)]){
              return ('Different')
            }else{
              return('Similar')
            }
          }
        }else{
          return('ERROR')
        }
      }else{
        return('Unshared')
      }
    }
  }))

  levelOrder<-factor(c('Identical','Shared','Different','Unshared'))
  a$shared<-factor(a$shared,levels=levelOrder)
  switch(typeR,
         'vec'= return(a$shared),
         'df' = return(a))
}


simpleVecSpe<-function(namedVec,breadthExp){
  namedVec<-sort(namedVec,decreasing=TRUE)
  namedVec<-breadthExp[names(namedVec)]==1
  newVec<-cumsum(namedVec)
  names(newVec)<-NULL
  return(newVec)
}

cleanupVec<-function(aVec,commonID){
  aVec<-aVec[commonID]
  aVec<-aVec[!is.na(aVec)]
  return(aVec)
}

