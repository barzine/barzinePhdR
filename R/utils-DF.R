# Headers manipulation -------------------------------

#' Prefix the names of the columns of a data.frame with a tag
#'
#' @param DF a data.frame
#' @param prefix a character vector; the tag to be added at the beginning
#'               of each column of the data.frame
#' @param sep the separator between the tag and the original content; default='.'
#'
#' @return a data.frame with the exact same content
#'         and where the columns have been prefixed with a given tag
#' @export
#'
prefix<-function(DF,prefix,sep='.'){
  if(missing(DF))
    stop('Missing argument: DF')
  if(missing(prefix))
    stop('Missing argument: prefix')
  if(!is.data.frame(DF))
    stop('object provided as DF is not from a compatible class')
  colnames(DF)<-sapply(colnames(DF), function(x) paste(prefix,x,sep=sep))
  return(DF)
}


#' Prefix character vector with a given tag
#'
#' @param vec a character vector
#' @param prefix the tag to be added at the beginning of each element of vec
#' @param sep the separator between the tag and the original content; default='.'
#'
#' @return a character vector
#' @export
#'
prefixv<-function(vec,prefix,sep='.'){
  if(missing(vec))
    stop('Missing argument: vec')
  if(missing(prefix))
    stop('Missing argument: prefix')
  if(!is.vector(vec))
    stop('object provided as vec is not from a compatible class')

  nvec<-sapply(vec,function(x) paste(prefix,x,sep=sep))
  return(nvec)
}

#' Tag each column name of a data.frame with a given string of characters
#'
#' @param DF a data.frame
#' @param suf the tag to be added at the end of each column name
#' @param sep the separator between the tag and the original content; default='.'
#'
#' @return a data.frame
#' @export
#'
suffix<-function(DF,suf,sep='.'){
  if(missing(DF))
    stop('Missing argument: DF')
  if(missing(suf))
    stop('Missing argument: suf')
  if(!is.data.frame(DF))
    stop('object provided as DF is not from a compatible class')

  colnames(DF)<-sapply(colnames(DF), function(x) paste(x,suf,sep=sep))
  return(DF)
}

# Conversion -----------------

#' Convert all the content of a given data.frame to characters
#'
#' @param DF a data.frame
#'
#' @return a character data.frame
#' @export
#'
factorDF2string<-function(DF){
  res <- data.frame(parallel::mclapply(DF, as.character), stringsAsFactors=FALSE)
  rownames(res)<-rownames(DF)
  return(res)
}


#' Convert a logical data.frame to a numeric data.frame (TRUE -> 1; FALSE -> 0)
#'
#' @param DF a boolean (logical) data.frame
#'
#' @return a numeric data.frame
#' @export
#'
boolDF2numeric<-function(DF){
  #res<-list()
  #res<-data.frame(lapply(colnames(DF),function(x) res[[x]]<-as.integer(DF[,x])))
  #rownames(res)<-rownames(DF)
  #colnames(res)<-colnames(DF)
  return(DF*1)
}


#' Change all factor columns in a data.frame into character column
#'
#' @param DF a data.frame (with a factor column)
#'
#' @return a data.frame
#' @export
#'
noMoreFactor<-function(DF){
  DF.rownames<-rownames(DF)
  DF<-data.frame(lapply(DF,function(x){
    if(is.factor(x))
      x<-as.character(x)
    return(x)
  }))
  rownames(DF)<-DF.rownames
  return(DF)

}


# Maths -------------------------------

#' Take two data.frames to create one and can log2 the content at once
#' while changing the -Inf to NA
#'
#' @param x a numeric data.frame
#' @param y a numeric data.frame
#' @param logT logical; default: FALSE
#'
#' @return one numeric data.frame
#' @export
#'
mkMatrix=function(x,y,logT=FALSE){
  if(logT){
    message("Logarithm: -Inf substitute by NA")
    return(as.matrix(apply(cbind(x,y),c(1,2),log2.na)))
  }else{return(as.matrix(cbind(x,y)))}
}

#' Replace one value by another and remove all NA rows
#' @description First intention meant to change -Inf to NA
#'              but should fonction with other changes as well
#'              as long as the values to be exchanged or from equivalent classes.
#'
#' @param DF a data.frame
#' @param motif the value that should be changed. Default: -Inf
#' @param token the new value. Default NA
#'
#' @return a data.frame
#' @export
#'
clean.infinite<-function(DF,motif=-Inf,token=NA){
  DF[DF==-Inf]<-token
  tmp<-stats::na.omit(DF)
  return(tmp)
}

## Randomisation

#' Randomise the content of each column of a data.frame
#'
#' @param DF a data.frame
#' @param seed to assure reproducibility. Default:213
#'
#' @return a randomised data.frame
#' @export
#'
dataset.randomize<-function(DF,seed=213){
  set.seed(seed)
  DF1<-as.data.frame(apply(DF,2,sample))
  row.names(DF1)<-row.names(DF)
  return(DF1)
}

#' Round numbers to the specified number of significant digits to keep
#'
#' @param DF numeric data.frame
#' @param i integer; indicates the number of significant digits to keep
#'
#' @return numeric data.frame
#' @export
#'
signifDF<-function(DF,i){
  DF.rownames<-rownames(DF)
  DF<-data.frame(lapply(DF,function(x){
    if(is.integer(x))
      return(x)
    if(is.numeric(x))
      x<-signif(x,digits=i)
    return(x)
  }))
  rownames(DF)<-DF.rownames
  return(DF)
}

## Stats ------------------------------

#' Arithmetical mean applied to a data.frame rows
#'
#' @param DF numeric data.frame
#'
#' @return the mean of each row (named vector)
#' @export
#'
RowMeans<-function(DF){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop("RowMeans need a numeric input"))
  res<-setNames(apply(DF,1,sum),rownames(DF))
  return(res)
}


#' Arithmetical median applied to a data.frame rows
#'
#' @param DF numeric data.frame
#'
#' @return the median of each row (named vector)
#' @export
#'
RowMedians<-function(DF){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop("RowMedians need a numeric input"))
  res<-setNames(apply(DF,1,stats::median),rownames(DF))
  return(res)
}

#' Median absolute deviation applied to a data.frame rows
#'
#' @param DF numeric data.frame
#' @param constant scale factor
#' @param ... parameters that can be passed to stats::mad
#'
#' @return the m.a.d of each row (named vector)
#' @export
#'
RowMad<-function(DF,constant,...){
  if(!missing(constant)){
    mad<-function(...){
      stats::mad(...,constant=constant)
    }
    return(mad)
  }
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop('RowMad need a numeric input'))
  corMAD<-setNames(apply(DF,1,mad),rownames(DF))
  return(corMAD)

}

#' Mode applied to a data.frame rows
#'
#' @param DF a numeric dataframe
#' @param type 'mode' for the mode; anything else gives the frequence of the mode
#'
#' @return the mode (or its frequence) of each row (named vector)
#' @export
#'
RowMode<-function(DF,type='mode'){
  return(setNames(apply(DF,1,function (x) {stat.mode(x,type)} ),rownames(DF)))
}


#' Standard deviation applied to a data.frame rows
#'
#' @param DF a numeric data.frame
#'
#' @return standard deviation of each row (named vector)
#' @export
#'
rowSD<-function(DF){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop('rowSD needs a numeric input'))
  corSD<-setNames(apply(DF,1,stats::sd),rownames(DF))
  return(corSD)
}

#' Variance applied to a data.frame rows
#'
#' @param DF a numeric data.frame
#'
#' @return variance of each row (named vector)
#' @export
#'
rowVar<-function(DF){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop('rowVar needs a numeric input'))
  corVar<-setNames(apply(DF,1,stats::var),rownames(DF))
  return(corVar)
}

#' Covariance applied to a data.frame rows
#' @description cv = sd/mean
#'
#' @param DF a numeric data.frame
#' @param digits integer indicating the number of significant digits
#'
#' @return the covariance of each row (named vector)
#' @export
#'
rowCV<-function(DF,digits=1){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop('rowCV needs a numeric input'))
  resCV<-setNames(signif(rowSD(DF)/rowMeans(DF),digits),rownames(DF))
  return(resCV)
}

## Covariance applied to a data.frame rows
## @description cv = sd/mean
##
## @param DF a numeric data.frame
##
## @return the covariance of each row (named vector)
## @export
##
#RowCV<-function(DF){
#  DF<-as.matrix(DF)
#  try(if(!is.numeric(DF)) stop('rowCV needs a numeric input'))
#  resCV<-setNames(signif(rowSD(DF)/rowMeans(DF)),rownames(DF))
#  return(resCV)
#}


#' Gives the maximum value of each row
#'
#' @param DF data.frame or matrix
#' @param type string to indicate how the result should be formatted
#'             "val" (default) returns the maximum values,
#'             "ind" returns the indices (as a named vector) of the maximum values,
#'             "all" returns a data.frame with the maximum value for each row and the indices.
#'
#' @return the maximum values and (or) their indices
#' @export
#'
rowMax<-function(DF,type='val'){
  switch(type,
         "val"={
           setNames(sapply(1:nrow(DF), function(x){
             max(DF[x,])
           }
           ),rownames(DF))
         },
         "ind"={
           sapply(1:nrow(DF), function(x){
             which.max(DF[x,])
           }
           )
         },
         "all"={
           tmp<-sapply(1:nrow(DF), function(x){ which.max(DF[x,]) })
           res<-as.data.frame(cbind(sapply(1:nrow(DF), function(x){ max(DF[x,])}),
                                    tmp,
                                    names(tmp)
           ),stringsAsFactors=FALSE)

           colnames(res)<-c("max",'indice','colname')
           rownames(res)<-rownames(DF)
           res[,1]<-as.numeric(res[,1])
           res[,2]<-as.integer(res[,2])
           res
         }
  )
}

#' Gives the minimum value of each row
#'
#' @param DF data.frame or matrix
#' @param type string to indicate how the result should be formatted
#'             "val" (default) returns the minimum values,
#'             "ind" returns the indices (as a named vector) of the minimum values,
#'             "all" returns a data.frame with the minimum value for each row and the indices.
#'
#' @return the minimum values and (or) their indices
#' @export
#'
rowMin<-function(DF,type='val'){
  switch(type,
         "val"={
           setNames(sapply(1:nrow(DF), function(x){
             min(DF[x,])
           }
           ),rownames(DF))
         },
         "ind"={
           sapply(1:nrow(DF), function(x){
             which.min(DF[x,])
           }
           )
         },
         "all"={
           tmp<-sapply(1:nrow(DF), function(x){ which.min(DF[x,]) })
           res<-as.data.frame(cbind(sapply(1:nrow(DF), function(x){ min(DF[x,])}),
                                    tmp,
                                    names(tmp)
           ),stringsAsFactors=FALSE)

           colnames(res)<-c("min",'indice','colname')
           rownames(res)<-rownames(DF)
           res[,1]<-as.numeric(res[,1])
           res[,2]<-as.integer(res[,2])
           res
         }
  )
}


#' Gives the maximum value of each column
#'
#' @param DF data.frame or matrix
#' @param type string to indicate how the result should be formatted
#'             "val" (default) returns the maximum values,
#'             "ind" returns the indices (as a named vector) of the maximum values,
#'             "all" returns a data.frame with the maximum value for each column and the indices.
#'
#' @return the maximum values and (or) their indices
#' @export
#'
colMax<-function(DF,type='val'){
  switch(type,
         "val"={
           setNames(sapply(1:ncol(DF), function(x){
             max(c(DF[,x]))
           }
           ),colnames(DF))
         },
         "ind"={
           tmp<-sapply(1:ncol(DF), function(x){
             which.max(c(DF[,x]))
           })
           return(setNames(tmp,rownames(DF)[tmp]))
         },
         "all"={
           tmp<-sapply(1:ncol(DF), function(x){
             which.max(c(DF[,x]))
           })
           res<-as.data.frame(cbind("max"=sapply(1:ncol(DF), function(x){ max(c(DF[,x]))}),
                                    "indice"=tmp,
                                    "rowname"=rownames(DF)[tmp]
           ),stringsAsFactors=FALSE)
           rownames(res)<-colnames(DF)
           class(res[,1])<-class(DF[,1])
           res[,2]<-as.integer(res[,2])
           res
         }
  )
}


#' Gives the minimum value of each column
#'
#' @param DF data.frame or matrix
#' @param type string to indicate how the result should be formatted
#'             "val" (default) returns the minimum values,
#'             "ind" returns the indices (as a named vector) of the minimum values,
#'             "all" returns a data.frame with the minimum value for each column and the indices.
#'
#' @return the minimum values and (or) their indices
#' @export
#'
colMin<-function(DF,type='val'){
  switch(type,
         "val"={
           setNames(sapply(1:ncol(DF), function(x){
             min(c(DF[,x]))
           }
           ),colnames(DF))
         },
         "ind"={
           tmp<-sapply(1:ncol(DF), function(x){
             which.min(c(DF[,x]))
           })
           return(setNames(tmp,rownames(DF)[tmp]))
         },
         "all"={
           tmp<-sapply(1:ncol(DF), function(x){
             which.min(c(DF[,x]))
           })
           res<-as.data.frame(cbind("min"=sapply(1:ncol(DF), function(x){ min(c(DF[,x]))}),
                                    "indice"=tmp,
                                    "rowname"=rownames(DF)[tmp]
           ),stringsAsFactors=FALSE)
           rownames(res)<-colnames(DF)
           class(res[,1])<-class(DF[,1])
           res[,2]<-as.integer(res[,2])
           res
         }
  )
}




#' Apply an Hampel's test on a data.frame
#'
#' @param DF a numeric data.frame
#' @param bool logical; default=TRUE. Whether the result should be given as a
#' @param threshold numeric. Default: 5.2 (expects a normal distribution)
#' @param ... additional parameters for stats::mad
#'
#' @return a result data.frame. If logical: TRUE means that the value passes the Hampel test for the given threshold.
#'                              If numeric, values failing the test are set to 0
#' @export
#'
hampel.test<-function(DF,bool=TRUE,threshold=5.2,...){
  DF<-as.matrix(DF)
  try(if(!is.numeric(DF)) stop('hampel.test need a numeric input'))
  if(bool){
    resDF<-setNames(lapply(rownames(DF),function(x){
      medx<-stats::median(unlist(DF[x,]))
      madx<-stats::mad(unlist(DF[x,],...))
      if(madx>0){
        return(abs(DF[x,]-medx)>t*madx)
      }else{
        return(0)
      }
    }),rownames(DF))
    resDF<-t(as.data.frame(resDF))
  }else{
    resDF<-setNames(lapply(rownames(DF),function(x){
      medx<-stats::median(unlist(DF[x,]))
      madx<-stats::mad(unlist(DF[x,],...))
      if(madx>0){
        return(abs(DF[x,]-medx)/madx)
      }else{
        return(0)
      }
    }),rownames(DF))
    resDF<-t(as.data.frame(resDF))
  }
  return(as.data.frame(resDF))
}

# Other ---------

#' Reverse a data.frame (or matrix) and change \code{-Inf} in \code{NA} or something else
#'
#' @param x numeric data.frame or matrix
#' @param replacement numeric, default NA.
#'
#' @return numeric data.frame
#' @export
#'
revCustom<-function(x,replacement=NA){
  x<-1/x
  x[x==Inf]<-replacement
  x[x==-Inf]<-replacement
  if(is.na(replacement)) x[is.na(x)]<-replacement
  return(x)
}

#' Reverse a symetric data.frame (or matrix). Specifically designed for symetric distance matrix.
#' Diagonal is identity, i.e. the distance on the diagonal is 0
#' (each object is the least distant to itself)
#'
#' @param x numeric data.frame or matrix
#' @param replacement numeric. Default 1 (to replace -Inf due to division by 0)
#'
#'
#' @return numeric data.frame
#' @export
#'
revCustomMatrixDist<-function(x,replacement=1){
  x<-1/as.matrix(x)
  diag(x)<-0
  x[x==Inf]<-replacement
  x[x==-Inf]<-replacement
  if(is.na(replacement)) x[is.na(x)]<-replacement
  return(x)
}



