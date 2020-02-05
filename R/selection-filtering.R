# Filtering on selection ------------------------------------------------------------
#' Selection on a data.frame to keep or filter out specified columns only
#'
#' @param DF data.frame to be filtered
#' @param IN column names to keep
#' @param OUT column names to filter out
#' @param ... other paramaters (for base::grep)
#' @param type pick the object to return;
#'             'DF' sends back the filtered dataframe
#'             'names' sends back the names (with the correct order) of the columns that are kept
#'
#' @return a dataframe or an ordered vector of the kept column names
#' @export
#'
colSelect<-function(DF,IN,OUT,...,type='DF'){

  colNis<-function(DF,pattern,...){
    return(colnames(DF)[grep(pattern,colnames(DF),...)])
  }

  colNisnot<-function(DF,pattern,...){
    return(colnames(DF)[-grep(pattern,colnames(DF),...)])
  }

  index<-rep(FALSE,ncol(DF))
  res<-unique(base::Reduce(c,Map(function(x) {grep(x,colnames(DF),...)} ,IN)))
  index[res]<-TRUE
  toRemove<-base::Reduce(c,Map(function(x) {grep(x,colnames(DF),...)},OUT))
  index[toRemove]<-FALSE
  switch(type,
         'names'={return(colnames(DF)[index])},
         'DF'   ={   new.DF<-DF[,index]
         return(new.DF)
         },
         #  'index'= #default (for the switch,function default: DF)
         {return(index)}
  )

}

#' Selection on a data.frame to keep the specified columns only
#'
#' @param DF data.frame to filter
#' @param IN names of the columns to be kept
#' @param ... other parameters (for base::grep)
#'
#' @return a filtered dataframe
#' @export
#'
colSelect2<-function(DF,IN,...){
  res<-base::Reduce(cbind, lapply(IN,function(x) DF[,grep(x,colnames(DF),...)]))
  res<-res[,unique(colnames(res))]
}

#' Keep only the specified rows of a dataframe
#'
#' @param DF data.frame to be filtered
#' @param filter vector of the index names of the rows to keep
#'
#' @return data.frame where the rows have been specified
#' @export
#'
filtre.dataset<-function(DF,filter){
  return(subDF<-DF[rownames(DF) %in% filter,]) #assignation eliminates the NA (?)
}


# Filtering based on exclusion -------------------------------
#' Create a data.frame for which specific rows have been excluded
#'
#' @param DF data.frame
#' @param rowlist names of the rows to be excluded
#'
#' @return filtered data.frame
#' @export
#'
excludeRows<-function(DF,rowlist){
  DF[!(rownames(DF) %in% rowlist),]
}

# Filtering based on numeric threshold ---------------------------------------
#' Test cell values in a data.frame for equality or higher value than
#' a specified cutoff
#'
#' @param DF numeric data.frame
#' @param CUT expression cutoff for the test
#'
#' @return a boolean data.frame where
#'         TRUE means that the value in the original data.frame is higher than CUT.
#' @export
#'
cutoff<-function(DF,CUT){
  apply(DF,2,function(x) x>=CUT)
}


#' Create a data.frame for which the values under a given threshold are set to 0
#'
#' @param DF numeric data.frame
#' @param CUT The values have to be higher than CUT to be kept, otherwise value is set to 0
#' @param clean if TRUE, the rows where no value are above cut are filtere out
#'
#' @return a data.frame which can have a smaller number of rows if clean is TRUE
#' @export
#'
df.cutoff<-function(DF,CUT,clean=TRUE){
  df1<-as.data.frame(apply(DF,c(1,2),function (x){
    if (x>=CUT){
      return(x)
    }else{
      return(0)}}))
  colnames(df1)<-colnames(DF)
  rownames(df1)<-rownames(DF)
  if (clean) {
    df1<-df1[rowSums(df1)!= 0, ]
  }
  return(DF[rownames(df1),])
}

#' Create a data.frame for which the values under a given threshold are set to 0
#'
#' @param DF numeric data.frame
#' @param cut the threshold
#' @param STRIP boolean; remove the rows filled with 0 when TRUE
#'
#' @return a numeric data.frame where all (value > cut|| value == 0)
#' @export
#'
df.cutoff2<-function(DF,cut,STRIP=FALSE){
  DFbool<-cutoff(DF,cut)
  DFbool<-boolDF2numeric(DFbool)
  res<-DFbool*DF
  if(STRIP){
    res<-stats::na.omit(res[rowSums(res)==0,])
  }
  return(res)
}

#' Strip the rows of a data.frame that are failling the test.
#'
#' @param DF data.frame to strip
#' @param test  'gt' (greater than), 'ge' (greater than or equal),
#'              'lt' (lesser than), 'le' (lesser than or equal)
#' @param threshold numeric value for the test
#' @param verbose when TRUE, shows the number of rows and columns before and after stripping
#'
#' @return a stripped data.frame
#' @export
#'
strip<-function(DF,test='gt',threshold=0,verbose=FALSE){
  if(verbose) {
    print(paste("Dimension of the dataset (before strip)",nrow(DF),"-",ncol(DF)))
  }
  switch(test,
         'gt'={
           DF1<-DF[rowSums(DF > threshold) > 0,]
         },
         'ge'={
           DF1<-DF[rowSums(DF %>=% threshold) >0,]
         },
         'lt'={
           DF1<-DF[rowSums(DF < threshold) > 0,]
         },
         'le'={
           DF1<-DF[rowSums(DF %<=% threshold) >0,]}
  )
  if(verbose) {
    print(paste("Dimension of the dataset (after strip)",nrow(DF1),"-",ncol(DF1),"(",test,"-",threshold,")"))
  }
  return(DF1)
}

#' Filter the content of a numeric data.frame based on a threshold
#'
#' @param DF the numeric data.frame
#' @param val the value used as a threshold
#' @param margin 1: for rows, 2 for columns, c(1:2) apply to both column and row.
#' @param case 'gt' for strictly greater than;
#'             'ge' for greater than or equal to;
#'             'lt' for strictly lesser than;
#'             'le' for lesser than or equal to;
#'             'eq' for equal to
#'
#' @return Data.frame which exclusively contains rows that comply at least once to the given criterion
#' @export
#'
df.Strip<-function(DF,val=0,margin=1,case='gt'){
  if(case=='gt')
    return(na.omit(DF[apply(DF,margin,sum)>val,]))

  if(case=='lt')
    return(na.omit(DF[apply(DF,margin,sum)<val,]))

  if(case=='ge')
    return(na.omit(DF[apply(DF,margin,sum)>=val,]))

  if(case=='le')
    return(na.omit(DF[apply(DF,margin,sum)<=val,]))

  if(case=='eq')
    return(na.omit(DF[apply(DF,margin,sum)==val,]))

  print("Error: not a correct case. Only possible cases: 'gt','ge','lt','le' or 'eq'")

}


#' Compare cell values to a threshold and fill cells
#' that failed the test with NA
#'
#' @param DF a numeric data.frame
#' @param val thereshold for the arithmetical comparison
#' @param case "ge","gt","le","lt","eq" (for df.Strip)
#'
#' @return a numeric data.frame
#' @export
#'
df.StripLine<-function(DF, val=0, case="gt"){
  df1<-df.Strip(DF, val=val, margin=c(1,2),case=case)
  df1[df1 == 0] <- NA
  df1<-stats::na.omit(df1[stats::complete.cases(df1),])
  return(df1)
}


#' Strip rows which don't have any value above the given threshold
#'
#' @param DF numeric data.frame
#' @param threshold threshold of detection
#'
#' @return a numeric data.frame which can have a smaller number of rows
#'         than the input data.frame
#' @export
#'
crudeStrip<-function(DF,threshold){
  DF<-apply(DF,c(1,2),function(x) {
    if(x<threshold){
      return(0)
    }else{
      return(x)
    }})
  DF<-strip(DF)
  return(as.data.frame(DF))
}

# Other filtering methods -------------------

#' Allow to select one tissue only from a data.frame
#'
#' @param DF a data.frame where the tissue are annotated as a part of the column names
#' @param Tissue character string
#'
#' @return a data.frame
#' @export
#'
createDFhXg<-function(DF,Tissue){
  newDF<-DF[,grep(Tissue,colnames(DF))]
  colnames(newDF)<-cleanNames(colnames(newDF),Tissue)
  return(newDF)
}

# Addition ----------------------------------
#' Create a data.frame with new rows or columns for the inputed data.frame
#' that are filled with a constant
#'
#' @param DF data.frame
#' @param NameList the header of the new rows or columns
#' @param value what the new cells in the data.frame should be filled with
#' @param type 'row' to create rows otherwise creates columns
#' @param check boolean
#' @param verbose print the different steps
#'
#' @return a completed data.frame (with the content of Namelist) and filled with the given value
#' @export
#'
df.fill<-function(DF,NameList,value=0,type='row',check=TRUE,verbose=TRUE){
  if(check){
    if (type=='row'){
      NameList<-base::setdiff(NameList,rownames(DF))
      print('Duplicates in new row names removed')
      if(length(NameList)==0) {
        if(verbose) print('Nothing to do')
        return(DF)
      }
    }else{
      NameList<-base::setdiff(NameList,colnames(DF))
      print('Duplicates in new column names removed')
      if(verbose) print('Nothing to do')
      return(DF)
    }
  }
  if(type=='row'){
    tDF<-t(as.matrix(DF))
    addDF<-base::list()
    addDF<-data.frame(lapply(NameList,function(x) addDF[[x]]<-rep(value,nrow(tDF))))
    rownames(addDF)<-rownames(tDF)
    colnames(addDF)<-NameList
    taddDF<-t(as.matrix(addDF))
    return(rbind.data.frame(DF,taddDF))
  }else{#every other choice creates columns
    addDF<-base::list()
    addDF<-data.frame(lapply(NameList,function(x) addDF[[x]]<-rep(value,nrow(DF))))
    rownames(addDF)<-rownames(DF)
    colnames(addDF)<-NameList
    return(cbind.data.frame(DF,addDF))
  }
}

