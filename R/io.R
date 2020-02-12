# Paths related ---------------------------------------

#if(!exists('dir.exists')) dir.exists <-function(paths){
#  file.exists(paste0(paths,'/'))
#}

#' Check if path is available for write
#' @param p path that has to be check for writing permission
#' @return adequate path
usablePath <- function(p) {
    try(if(!dir.exists(p))
        stop(paste(p, 'not available as writing directory')))

    if(substr(p,nchar(p),nchar(p))!='/')
        p=paste(p,'/',sep='')
    return(p)
}

# Input -----------------------------------------------------

#' Basically a read.table but with predetermnined parameters
#' @param filename Name of the file of the data.frame to read
#' @param path Path of the file if not included in the name or not the working directory
#' @param row.names indicates the row that contains the column names
#' @param ... any other parameter supported by utils::read.table
#' @return an adequately formatted data.frame
#' @export
Read.table <- function(filename,...,path, row.names=1){
  if (missing(path) &
      !as.logical(length(grep('/',filename))) &
      !exists('read.PATH', envir=globalenv()))
    path='./'
  if(as.logical(length(grep('/',filename))))
    path=''
  filename=paste0(path,filename)
  return(utils::read.table(filename,stringsAsFactors=FALSE,
                           header=TRUE,sep='\t',row.names=row.names,...))
}


#' Create a vector from a file where the data is formatted as a list
#' with each element on a newline
#'
#' @param filename path to the file
#'
#' @return a vector
#' @export
#'
listFromFile<-function(filename){
  l<-read.table(filename,header=FALSE, stringsAsFactors=FALSE)
  l<-l[,1]
  l<-unique(l)
  return(l)
}

#' Create a vector from a file downloaded from Biomart
#' where the data is formatted as a list
#' each newline comprises a gene id and transcript id.
#'
#' @param filename path to the file
#' @param out boolean. Default:TRUE. Whether the genes ID should be returned
#'
#' @return a vector of gene ID or a data.frame with the genes and transcripts ID
#' @export
#'
listFromFileBiomart<-function(filename,out=TRUE){
  l<-read.table(filename,header=TRUE, stringsAsFactors=FALSE,sep='\t')
  if(out){
    l<-l[,1]
    l<-unique(l)
    return(l)
  }else{
    if(dim(l)[2]==2){
      colnames(l)<-c('Gene.ID','Transcript.ID')
    }

  }
  return(unique(l))
}

#' Retrieve all the genes ID from a file
#'
#' @param filename path to the file
#' @param h character string. Header of the files. Default: "Gene.ID"
#' @return a vector or a named vector if there are headers
#' @export
listgenes<-function(filename,h="Gene.ID"){
  l<-read.table(filename,header=TRUE, stringsAsFactors=FALSE)
  l<-l[,h]
  l<-unique(l)
  return(l)
}

#' Read a file (one or two columns) and creates a vector (named if two-column file)
#'
#' @param filename Name of the file
#' @param header boolean, TRUE if the file contains the header of the data
#' @param stringsAsFactors default FALSE
#' @param trunc default FALSE
#' @param ... other parameters can be passed to utils::read.table
#'
#' @return a vector or a named vector if there are headers
#' @export
#'
read.list<-function(filename, header=TRUE, stringsAsFactors=FALSE, trunc=FALSE,...){
  tmp<-utils::read.table(filename,header=header, stringsAsFactors=FALSE,...)
  vec<-tmp[,1]
  names(vec)<-rownames(tmp)
  if(trunc){
    vec<-vec[vec>0]
  }
  return(vec)
}

#' Load function that returns the object(s) instead of attaching it to the global namespace
#' stolen from K. Rudolph and M. Schubert in the `ebits` repo
#' @param filename filename to be loaded in the global environment
#' @return the object content which needs to be stored for later use
#' @export
Load <- function(filename) {
  lfc = function(fpath) {
    env = new.env()
    fdir = dirname(fpath)
    fid = strsplit(basename(fpath), "\\$")[[1]]
    fname = fid[1]
    subsets = fid[-1]

    base::load(file.path(fdir,fname), env)
    contents = as.list(env)
    if (length(contents)==1)
      contents[[1]]
    else
      contents
  }
  if (length(filename) > 1)
    lapply(filename, lfc)
  else
    lfc(filename)
}

#' Load in the working environment the object matching the name given as input
#' (will look for DIR.PATH in the environment)
#' can only load RData files
#' @param objectname Name of the .RData file (without the extension)
#' for which an object has to be created in the environment
#' @param path Path of the object to be loaded
#' @param env environment where the object should be created
#' Default is global environment
#' @return Output directly the object inside the specified environment
#' @export
qLoad<-function(objectname, path, env=globalenv()){
  if (missing(path) & !as.logical(length(grep('/',objectname))) & !exists('DIR.PATH', envir=globalenv())){
    path <- './'
  }else{
    if(missing(path) & exists('DIR.PATH', envir=globalenv())){
      path <- globalenv()$DIR.PATH
    }
  }
  return(base::load(paste0(path,objectname,'.RData'),envir=env))
}

# Output --------------------------------------------------------------
#' Takes an object and saves it in a file with the same name.
#' Has an alternative where the name can be changed: saveWithName
#' @param obj the object to be saved
#' @param filename the filename can be specified. If it is not, the name object will be used
#' @param path the path where the object should be saved
#' @param overwrite if a file exists at the specified, if it should be overwrite with the new object
#' @param extension preferable to not change. By default: '.Rdata'
#' @param show.path Either if the path should be shown in the message. Default: `FALSE`
#' @return message that indicates if the object has been saved (then with its name)
#' @export
saveToFile<-function(obj,filename,path,overwrite,extension='.RData', show.path=FALSE){
  if(missing(filename))
    filename=substitute(obj)
  if(missing(path)){
    if(!exists('write.PATH', envir=globalenv())){
      path <- './'
    }else{
      path <-globalenv()$write.PATH
    }
  }
  path <- usablePath(path)

  if(missing(overwrite)){
    if(exists('overwrite', envir=globalenv())){
      overwrite<-globalenv()$overwrite
    }else{
      overwrite<-FALSE
    }
  }

  if(length(grep('.RData',filename)>0)&&grep('.RData',filename)==1) extension=''

  out <- paste0(path,filename,extension)

  obj <- stats::setNames(list(obj),deparse(substitute(obj)))

  if(!file.exists(out)){
    save(list=names(obj),file=out,envir=list2env(obj))
    message <- paste0(names(obj),' save as ',filename,extension)
  }else{
    if(!overwrite){
      message <- paste(paste0(filename,extension),"exists already at the specified path - nothing done")
    }else{
      message=paste('over writting',out)
      save(list=names(obj),file=out,envir=list2env(obj))
      message <- paste(message,paste(filename,'save as',filename,extension),sep='\n')
    }
  }
  if(show.path) message <-paste(message,'path:',path)
  return(message)
}


#' Allows changing the name of the object at the same time that it is saved to a file
#' @param ... object to be saved
#' @param file new name of the object which is saved with the same name
#' @return saved file
#' @export
saveWithName <- function(...,file){
  obj<-list(...)
  save(list=names(obj),file=file,envir=list2env(obj))
}


#' write.table applied to DF plus some general customisation
#' @param DF data.frame to be saved to a text file
#' @param filename Name of the file to be saved
#' @param path Path where the file should be saved
#' @param extension by default '.tsv'
#' @param overwrite if `TRUE` indicates that the file should be saved in any case,
#' even if it overwrites a previous file with the same name at the same path
#' @param show.path indicates if the path in the message should be printed too
#' @param rowNames The name of the first column in the file
#' @return message that indicates if the object has been saved to a (named) text file
#' @export
wtDF<-function(DF,filename,path,
               rowNames="Gene.ID.ens76",overwrite,
               extension='.tsv', show.path=FALSE){
  if(missing(overwrite)){
    if(exists('overwrite', envir=globalenv())){
      overwrite<-globalenv()$overwrite
    }else{
      overwrite<-FALSE
    }
  }
  if(exists("write.PATH", envir=globalenv())) {
    path<-usablePath(globalenv()$write.PATH)
  }else{
    path<-'./'
  }

  if(missing(filename)) filename<-substitute(DF)

  filename<-paste0(filename,extension)
  out<-paste0(path,filename)

  new.DF<-data.frame(rownames(DF),DF)
  colnames(new.DF)[1]<-rowNames


  if(!file.exists(out)){
    utils::write.table(new.DF,file=filename,row.names=FALSE,sep='\t')
    message <- paste0(substitute(DF),' save as ',filename)
  }else{
    if(!overwrite){
      message <- paste(filename,"exists already at the specified path- nothing done")
    }else{
      message=paste('over writting',filename)
      utils::write.table(new.DF,file=filename,row.names=FALSE,sep='\t')
      message <- paste(substitute(DF),'save as',filename)
    }
  }
  if(show.path) message<-paste(message,'- path:',path)
  return(message)
}


# Helper functions -------------------------------------

#' Convert all missing values to 0
#'
#' @param DF numeric data.frame or vector
#' @param missingVal value to change to 0
#'
#' @return a data.frame with no missing value
#' @export
#'
NAto0<-function(DF,missingVal=NA){
  DF[is.na(DF)]<-0
  return(DF)
}

#' Creates a two-column data.frame from a named vector where the column names are repeated in the first column
#'
#' @param filename Name of the file which comprises a named vector
#'
#' @return a data.frame of two columns where the first column is the repetition of the row.names
#' @export
#'
nl1d<-function(filename){
  DF<-utils::read.table(filename,header=FALSE, stringsAsFactors=FALSE)
  rownames(DF)<-DF[,1]
  colnames(DF)<-c('ID',strsplit(basename(filename),"\\.")[[1]][1])
  return(DF)
}







