#' Format a data.frame in long format and annotate columns
#'
#' @param DF a data.frame
#' @param RowNames weither the rownames should be included as new column
#' @param TRow Name of the column for the row.names. Default="Gene".
#' @param categoriesA vector comprising the different string that should be look into
#'                    the column names of the input data.frame
#' @param categoriesB vector comprising the different string that should be look into
#'                    the column names of the input data.frame
#' @param nameA Name of the new column to add for categoriesA in the result data.frame
#' @param nameB Name of the new column to add for categoriesB in the result data.frame
#' @param col1 the name of the column to parse for the categories.
#'            default: "variable" (default name given by reshape2::melt for expression data.frame )
#'
#' @return a long data.frame
#' @export
#'
lg.format<-function(DF,RowNames=TRUE,TRow="Gene",categoriesA,categoriesB,
                    nameA="condition",nameB="study",col1="variable"){

  if(RowNames){
    DF$new<-rownames(DF)
    colnames(DF)[ncol(DF)]<-TRow
  }

  newDF<-reshape2::melt(DF)

  if(!missing(categoriesA)){
    newDF$nameA<-NA
    newDF<-base::Reduce(rbind,parallel::mclapply(categoriesA,function(x){
      tmp<- newDF[grep(x,newDF[,col1]),]
      tmp$nameA<-x
      stats::na.omit(tmp)
    })
    )

    colnames(newDF)[grep('nameA',colnames(newDF))]<-nameA
  }

  if(!missing(categoriesB)){
    newDF$nameB<-NA
    newDF<-base::Reduce(rbind,parallel::mclapply(categoriesB,function(x){
      tmp<- newDF[grep(x,newDF[,col1]),]
      tmp$nameB<-x
      stats::na.omit(tmp)
    })
    )

    colnames(newDF)[grep('nameB',colnames(newDF))]<-nameB
  }
  return(newDF)
}

#' Create a data.frame from the elements of two vectors that have the same names
#'
#' @param a first vector
#' @param b second vector
#' @param colN Colnames for the new data.frame
#'
#' @return a data.frame with two colums
#' @export
#'
make.2ldf<-function(a,b,colN){
common<-intersect(names(a),names(b))
tmpa<-a[common]
tmpb<-b[common]

DF<-cbind(tmpa,tmpb)
if(!missing(colN)) colnames(DF)<-colN
return(DF)
}

#' Create a data.frame from the elements of three vectors that have the same names
#'
#' @param a first vector
#' @param b second vector
#' @param c third vector
#' @param colN Colnames for the new data.frame
#'
#' @return a data.frame with three columns
#' @export
#'
make.3ldf<-function(a,b,c,colN){
 common<-intersect(intersect(names(a),names(b)),names(c))
 tmpa<-a[common]
 tmpb<-b[common]
 tmpc<-c[common]

 if(!missing(colN)) DF<-cbind(tmpa,tmpb,tmpc)
 colnames(DF)<-colN
 return(DF)
}


#' Selection of the same conditions (columns) between two studies
#'
#' @param df1 first data.frame
#' @param df2 second data.frame
#' @param cond character string vector.
#'             A list can be given that allows to make a preselection on the kept conditions.
#' @param name1 character string. Allows to prefix the columns names of the first data.frame
#' @param name2 character string. Allows to prefix the columns names of the second data.frame
#'
#' @return one data.frame with the selected shared conditions of the two data.frames
#' @export
#'
cross.selectCond<-function(df1,df2,cond=NA,name1=NA,name2=NA){
  #to avoid futur errors the conditions are checked:
  if (!is.na(cond[1])){
    cond<-intersect(colnames(df1),cond)
    cond<-intersect(colnames(df2),cond)
  }else{
    cond=intersect(colnames(df1),colnames(df2))
  }
  if (is.na(name1)) 	name1=deparse(substitute(df1))
  if (is.na(name2)) 	name2=deparse(substitute(df2))
  return(Reduce(cbind,lapply(cond,function(x) {
    select=intersect(row.names(df1),row.names(df2))
    df<-cbind(df1[select,x],df2[select,x])
    rownames(df)<-select
    colnames(df)<-c(paste(name1,x,sep='.'),paste(name2,x,sep='.'))
    return(as.data.frame(df))
  }
  )))
}


#' Keep only the columns or rows or both that are common between the two data.frames
#'
#' @param df1 data.frame for the first study
#' @param df2 data.frame for the second study
#' @param margin character. Default: 'c'.
#'               Selection only on the columns "c",
#'               selection only on the rows "r",
#'               selection on both "a".
#' @return a list of two data.frames
#' @export
#'
cross.select<-function(df1,df2,margin='c'){
  if (!exists("sharedEnv")) barzinePhdR::initialise()


  #check that given margin is allowed
  if (!margin %in% c('c','r','a')) stop("Margin should be 'c', 'r' or 'a'")

  if (margin=='c' || margin=='a'){ #reduction on columns
    common<-intersect(colnames(df1),colnames(df2))
    if (length(common)>0){
      dfa<-df1[,common]
      dfb<-df2[,common]
    }else{#there is not any equivalent attribute between df1 and df2
      return(list(NA,NA))
    }

    if (margin=='a'){ #selection on rows
      common<-intersect(rownames(dfa),rownames(dfb))
      if (length(common)>0){
        dfa<-dfa[common,]
        dfb<-dfb[common,]
      }else{#there is not any equivalent attribute between df1 and df2
      }
      return(list(NA,NA))
    }
    return(list(dfa,dfb))
  }

  if (margin=='r'){ #selection on rows
    common<-intersect(rownames(df1),rownames(df2))
    if (length(common)>0){
      dfa<-df1[common,]
      dfb<-df2[common,]
    }else{#there is not any equivalent row between df1 and df2
      return(list(NA,NA))
    }

    return(list(dfa,dfb))
  }
  print('ERROR in cross.select()')
}

#' Allow to save a data.frame quicly.
#'
#' @param ... the data.frame to save
#' @param path character string. Path where to save the data.frame
#' @param file character string. Name of the file.
#'             If left empty, a time stamped is used instead.
#'
#' @return the name of the saved file
#' @export
#'
save_df<-function(...,path='../data/Robject', file=NA){
  if(is.na(file[1])){
    #the name of the file would be the date of the save
    file<-paste('save',format(Sys.time(), "save_%d-%b-%Y_%H-%M-%S"),sep="")

  }else{
    file<-paste(file,format(Sys.time(), "%d-%b-%Y_%H-%M-%S"),sep='_')
  }
  file<-paste(file,'.RData',sep="")
  save(...,file=paste(path,file,sep='/'),ascii = TRUE)
  print('Save done!')
  print(file)
}

#' Attach the current time to a string name
#'
#' @param string character string
#'
#' @return the string with the current time
#' @export
#'
stamped<-function(string){
  return(paste(string,format(Sys.time(), "save_%d-%b-%Y_%H-%M-%S"),sep=""))

}

#' Create a data.frame for a transcriptomic study from the different files outputted by older versions of irap
#'
#' @param path path to the directory that contains all the files for a given study
#'
#' @return a data.frame
#' @export
#'
dataset_irap<-function(path){

  files<-as.list(list.files(path=path, full.names=TRUE))
  df1<-lapply(files,extract.irap)
  #merge the different lists and reduce the dataframe
  df2<-Reduce(merge, df1)
  #copy the genes names to the rownames of the data.frame
  rownames(df2)<-df2$ID
  #delete the "ID" column
  df2<-df2[,!names(df2) %in% "ID"]

  return(df2)
}

#' Import data from files produced by an older version of irap.
#'
#' @param filename character string. Path/filename of a sample (?)
#' @return a data.frame
#' @export
#'
extract.irap<-function(filename){
  df.raw<-read.table(filename,header=TRUE, stringsAsFactors=FALSE)
  df<-df.raw[df.raw$FPKM_status=='OK',c('Gene','FPKM')]
  colnames(df)<-c('ID',strsplit(basename(filename),"\\.")[[1]][1])
  #cleanup of the duplicated
  df1<-aggregate(. ~ ID, data=df, FUN=sum)
  rownames(df1)<-df1$ID
  return(df1)
}


#' Wrapper around read.table with a few tweaks
#'
#' @param file path to the file
#' @param header boolean. Default: TRUE. Whether the file contains headers.
#' @param sep character. Default: "\ t"
#' @param stringsAsFactors boolean. Default=FALSE.
#' @param ... other arguments for read.table
#' @param select Value to keep (?)
#' @param rnames default:NA. vector of rownames to keep (?)
#' @param drops default:NA. vector of column names that should be removed from the data.frame
#'
#' @return a data.frame
#' @export
#'
creationDf.gxa<-function(file,header=TRUE,sep="\t",stringsAsFactors=FALSE,...,
                         select=NA,rnames=NA, drops=NA){
  df<-read.table(file=file,header=header,sep=sep,stringsAsFactors=stringsAsFactors,...)
  if(!is.na(rnames[1])){
    row.names(df)<-df[[rnames]]
  }
  if(!is.na(select[1])){
    df=df[select]
  }
  if(!is.na(drops[1])){
    df<-df[,!(names(df) %in% drops)]
  }

  if(row.names(df)[1]==1){
    print("Warning: The row names don't seem to be assigned. Check that the headers of the file match with the given argument rnames")
  }

  return(df)
}


#' Create a data.frame from the two data.frames after filtering the columns and rows.
#'
#' @param df1 data.frame for a first study
#' @param df2 data.frame for a second study
#' @param cond1 vector of strings of character. Names of the columns to keep for df1
#' @param cond2 vector of strings of character. Names of the columns to keep for df2.
#'              default: "cond2=cond1"
#' @param select vector of strings of character. Names of the rows to keep in df1 and df2.
#'               default: intersect(row.names(df1),row.names(df2)))
#'
#' @return a data.frame created from the junction of the df1 and df2
#'         after the filtering.
#' @export
#'
refGraphDF<-function(df1,df2,cond1,cond2=cond1,
                     select=intersect(row.names(df1),row.names(df2))){
  df<-cbind(df1[select,cond1],df2[select,cond2])
  rownames(df)<-select
  colnames(df)<-c(paste(deparse(substitute(df1)),cond1,sep='.'),
                  paste(deparse(substitute(df2)),cond2,sep='.'))
  return(as.data.frame(df))
}

#' Remove rows from data.frame (created by older version of irap)
#' that contain additionnal inforamtion
#'
#' @param df data.frame created from older version of irap
#'
#' @return data.frame that contains gene expression only
#' @export
#'
clean.HTSc<-function(df){
  df<-df[rownames(df) !="alignment_not_unique",]
  df<-df[rownames(df) !="not_aligned",]
  df<-df[rownames(df) !="too_low_aQual",]
  df<-df[rownames(df) !="ambiguous",]
  df<-df[rownames(df) !="no_feature",]
}


#' Allows to format older version of irap's output
#'
#' @param filename path to the file
#' @param reduce boolean. Default: TRUE. Whether the "ID" column should be removed (duplicated with the rownames)
#' @param clean boolean. Default:TRUE.
#'
#' @return data.frame of gene expression
#' @export
#'
extract.HTSc<-function(filename,reduce=FALSE,clean=TRUE){
  df<-read.table(filename,header=FALSE, stringsAsFactors=FALSE)
  colnames(df)<-c('ID',strsplit(basename(filename),"\\.")[[1]][1])
  rownames(df)<-df$ID

  if (clean) {
    #the last lines in the dataset are removed
    df<-clean.HTSc(df)
  }
  if(reduce){
    #the ID column is removed
    df$ID<-NULL
  }

  return(df)
}


#' Create a data.frame from a folder of htseq outputs
#'
#' @param pathFolder path to the directory
#' @param clean boolean. Default: TRUE. Whether the additional annotations should be removed from the result data.frame
#'
#' @return an expression (numeric) data.frame
#' @export
#'
dataset.HTSc<-function(pathFolder,clean=TRUE){
  files<-list.files(path=pathFolder, full.names=TRUE)
  df1<-lapply(files,nl1d)
  #merge the different lists and reduce the dataframe
  df2<-Reduce(function(x,y) {merge(x,y, by='ID')}, df1)
  df2<-df2[,!names(df2) %in% "ID"]
  rownames(df2)<-df1[[1]]$ID

  if (clean) {
    return(clean.HTSc(df2))
  }else{
    return(df2)
  }
}


#' Create a dataframe filled with a value with the same number of rows than a given df
#' and then bind the two together (new data.frame as new columns)
#'
#' @param df data.frame of reference
#' @param att vector of character strings
#' @param val vector of values to use to fill the new columns
#' @param ... additionnal arguments for base::data.frame
#'
#' @return a data.frame
#' @export
#'
multi.add.qual.att<-function(df,att,val,...){
  l<-lapply(val,times=dim(df)[1],rep)
  dftmp<-data.frame(Reduce(cbind,l),...)
  colnames(dftmp)<-att
  rownames(dftmp)<-rownames(df)
  return(cbind(df,dftmp))

}

#' Create a function that allows to add a qualitative attribute easily
#'
#' @param att attribute that is fixed in the new function
#'
#' @return a function
#' @export
#'
gen.add.qual.att<-function(att){
  function(df,val){
    df[dim(df)[2]+1]<-val
    colnames(df)[dim(df)[2]]<-att
    return(df)
  }
}

#' Add a (named) column to a data.frame, filled with the same value
#'
#' @param df data.frame
#' @param att Name of the new attribute (name of the new column)
#' @param val value with which fill the new column
#'
#' @return data.frame
#' @export
#'
add.qual.att<-function(df,att,val){
  df[dim(df)[2]+1]<-val
  colnames(df)[dim(df)[2]]<-att
  return(df)
  #note: other way to create the dataframe:
  #df<-date.frame(df,val);colnames(df)[dim(df)[2]]<-att
}


#' From a file outputted by an older version of irap,
#' create a data.frame and add new columns with a selected generic value
#'
#' @param filename character strings. Path of the file
#' @param att vector of attributes to add
#' @param value vector with the values with which to fill the new columns
#' @param ... other arguments handled by multi.add.qual.att
#'
#' @return a data.frame
#' @export
#'
single.factoring<-function(filename,att,value,...){
  predf<-extract.HTSc(filename,...)
  return(multi.add.qual.att(predf,att=att,val=value,...))
}



#' From a filename, extract a data.frame to which a number of columns
#' (filled with selected values )
#'
#' @param filename character string, path to the file
#' @param att a vector of character strings (names of the new columns)
#' @param value vector of values for the new columns
#' @param ... other parameters that can be handled by multi.add.qual.att
#'
#' @return a data.frame
#' @export
#'
single.factoring2<-function(filename,att,value,...){
  predf<-extract.HTSc(filename,...)
  tmp<-colnames(predf)[2]
  colnames(predf)[2]<-'counts'
  att1<-c('sample',att)
  val1<-c(tmp,value)
  return(multi.add.qual.att(predf,att=att1,val=val1,...))

}

#' From a list of filenames, return a list of data.frames after adding an attribute
#'
#' @param att1 an attribute
#' @param data filenames
#' @param ... other arguments for single.factoring
#'
#' @return a list of data.frames
#' @export
#'
factoring<-function(att1,data,...){
  listDF<-lapply(data,function(x) {
    file<-x[['filename']]
    val<-x[['value']]

    single.factoring(filename=file,att=att1,value=val,...)
  })
  return(list(listDF))
}

#' From a list of data.frame, add an attribute and collate (row binding) them
#'
#' @param att1 attribute to add to the data.frames
#' @param data list of data.frame
#' @param ... other arguments for single.factoring2
#'
#' @return a data.frame
#' @export
#'
factoring.rbind<-function(att1,data,...){
  listDF<-lapply(data,function(x) {
    file<-x[['filename']]
    val<-x[['value']]

    single.factoring2(filename=file,att=att1,value=val,...)
  })
  return(Reduce(rbind,listDF))
}






