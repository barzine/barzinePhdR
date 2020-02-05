#' Allows to load an annotation data.frame
#'
#' @param PATH path of the annotation file
#' @param standardise boolean; default: FALSE.
#'                    Whether the tissue names should be easier to matched
#'                    across datasets.
#'
#' @return a data.frame that contains the annotation data associated to some of the
#'        dataset used in the PhD project
#' @export
#'
cleanReadAnnot<-function(PATH,standardise=FALSE){
  if(tools::file_ext(PATH) =='tsv'){
    tmp<-utils::read.table(PATH, stringsAsFactors=FALSE,header=TRUE, sep="\t")
    rownames(tmp)<-tmp$sample
  }else{
    if(tools::file_ext(PATH)=='RData') tmp<-Load(PATH)
  }

  if(standardise){
    tmp$condition<-gsub('Cerebral.cortex','Cortex',tmp$condition)
    tmp$condition<-gsub('Endometrium','Uterus',tmp$condition)
    tmp$condition<-gsub('Salivary.gland','Salivary',tmp$condition)
    tmp$condition<-gsub('Spinal.cord','Spinalcord',tmp$condition)
    tmp$condition<-gsub('Gall.bladder','Gallbladder',tmp$condition)
  }
  return(tmp)
}


#' Wrapper that briefly present the datasets
#'
#' @param DF a numeric (expression) dataset
#' @param label.df string of character that allows to label the dataset in the output
#' @param feature string of characters (essentially either 'gene', 'mRNA', 'RNA' or 'protein')
#'
#'
#' @return a brief summary for the dataset
#' @export
#'
presentData<-function(DF,label.df,feature){
  cat(paste("###",label.df,'\n\n'))

  check<-nbGenesPerCond.unique(DF,
                               xlab='Tissue',
                               title=paste0("Distribution of unique ",feature," per tissue \n for ",label.df))
  plot_distribDF(DF,colourVal=barzinePhdData::TissueCol_21,verbose=FALSE)

  return(check)
}

#' Check that both datasets share identical dimensions
#'
#' @param DF1 first data.frame
#' @param DF2 second data.frame
#'
#' @return summary of the conformity of the two dataset or not
#' @export
#'
checkConformity<-function(DF1,DF2){
  if (!all(dim(DF1)==dim(DF2)))  print(dim(DF1)==dim(DF2))
  if (!all(rownames(DF1)==rownames(DF2))) print(rownames(DF1)==rownames(DF2))
  if (!all(colnames(DF1)==colnames(DF2))) print(colnames(DF1)==colnames(DF2))

  print(paste('Conformity check done for',substitute(DF1),'and',substitute(DF2)))
}


#' Load three data.frames from file paths and check that they have the same dimensions
#'
#' @param path1 file path to the first data.frame
#' @param path2 file path to the second data.frame
#' @param path3 file path to the third data.frame
#'
#' @return a list of three data.frames
#' @export
#'
load_and_check_3DF<-function(path1,path2,path3){

  #Load
  DF1=Load(path1)
  DF2=Load(path2)
  DF3=Load(path3)

  #warning if there are any NA:
  if(any(is.na(DF1))) warning(paste("There are missing values in (DF1)",path1))
  if(any(is.na(DF2))) warning(paste("There are missing values in (DF2)",path2))
  if(any(is.na(DF3))) warning(paste("There are missing values in (DF3)",path3))

  #check that the dimensions are the same
  checkConformity(DF1,DF2)
  checkConformity(DF1,DF3)
  checkConformity(DF2,DF3)
  return(list(DF1,DF2,DF3))
}


#' Rows that are unmatched in one dataset by any row in two other datasets.
#'
#' @param a first data.frame
#' @param b second data.frame
#' @param c third data.frame
#'
#' @return Names of the unmatched rows
#' @export
#'
unmatchedRowsAbyBC<-function(a,b,c){
  a<-rownames(a)
  b<-rownames(b)
  c<-rownames(c)
  return(
    setdiff(a,Union(Intersect(a,b,c),
                    intersect(a,b),
                    intersect(a,c)
    )
    )
  )
}

#' Rows that are unmatched in one dataset by any row in another dataset.
#'
#' @param a first data.frame
#' @param b second data.frame
#'
#' @return Names of the unmatched rows
#' @export
#'
unmatchedRowsAbyB<-function(a,b){
  a<-rownames(a)
  b<-rownames(b)
  return(setdiff(a,intersect(a,b)))
}


