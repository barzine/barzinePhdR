#' Allows to output a list directly in several objects (python like behaviour)
#' @description This function is interesting with functions
#' that handle several objects at the same time
#' The function was created by Gabor Grothendieck
#' who provided it graciously on the r-help list in June 2004.
#' https://stat.ethz.ch/pipermail/r-help/2004-June/053343.html
#' He originally named the function list (overriding then the built-in function).
#' However, I have preferred to rename it for avoiding confusion.
#'
#' @param x list of objects
#' @param ... other parameters
#' @param value value
#'
#' @return a list of object that can be outputed directly in a list of objects
#' @export
#'
"[<-.result" <- function(x,...,value) {
  args <- as.list(match.call())
  args <- args[-c(1:2,length(args))]
  length(value) <- length(args)
  for(i in seq(along=args)) {
    a <- args[[i]]
    if(!missing(a)) eval.parent(substitute(a <- v,list(a=a,v=value[[i]])))
  }
  x
}

# Maths -------------------------------

sharedEnv<-new.env(parent = emptyenv())
assign("epsilon",1e-10,envir = sharedEnv)

#' Allows to change the value of epsilon (needed to circumvent some R caveats)
#'
#' @param x New value for epsilon; default 1e-10
#'
#' @return Change the value of epsilon used by the other functions
#' @export
#'
#' @examples setEpsilon(1e-9)
setEpsilon<-function(x){
  assign("epsilon",x,sharedEnv)
  message(paste('Epsilon set to',sharedEnv$epsilon))
}


#' Displays the value of epsilon
#'
#' @return the value of epsilon
#' @export
#'
getEpsilon<-function(){
  message(paste('Epsilon is set to',sharedEnv$epsilon))
}

#' Test more carefully if x is equal or greater than y
#'
#' @param x first member of the comparison
#' @param y second member of the comparison
#'
#' @return boolean
#' @export
#'
#' @examples 8 %>=% 8.0
`%>=%` <- function(x, y){
  x + sharedEnv$epsilon > y
}

#' Test more carefully if x is equal or lesser than y
#'
#' @param x first member of the comparison
#' @param y second member of the comparison
#'
#' @return boolean
#' @export
#'
#' @examples 8 %<=% 8.00000001
`%<=%`<- function(x, y){
  x - sharedEnv$epsilon < y
}

#' Statistical mode value
#' @description Compute the sample mode
#'              Nearly copy-paste from an answer from Ken Williams http://stackoverflow.com/users/169947/ken-williams
#'              see http://stackoverflow.com/a/8189441/2116422
#'
#' @param x a numeric vector
#' @param type mode' for the mode; anything else gives the frequence of the mode
#'
#' @return the mode (single numeric value) or its frequence
#' @export
#'
stat.mode <- function(x,type='mode') {
  if(type=='mode'){
    unique_values <- unique(x)
    return(unique_values[which.max(tabulate(match(x, unique_values)))])
  }else{ #return the frequence of the mode
    return(max(tabulate(match(x, unique_values))))
  }
}

#' The correct answer -Inf for log2(0) is changed to NA
#'
#' @param x numeric value
#'
#' @return log2(x) or NA for log2(O)
#' @export
#'
log2.na=function(x) {
  if(log2(x)==-Inf){
    return(NA)
  }else{
    return(log2(x))
  }
}

# Set operations ------------------

#' Allow to intersect more than two sets at a time
#'
#' @param ... comma-separetated list of sets
#'
#' @return the intersection of all the sets
#' @export
#'
Intersect <- function(...) {
  base::Reduce(base::intersect, list(...))
}

#' Allow to intersect more than two sets at a time
#'
#' @param ... comma-separetated list of sets
#'
#' @return the union of all the sets
#' @export
#'
Union <- function(...) {
  base::Reduce(base::union, list(...))
}

#have to be tested!
#Diff<-function(...){Reduce(base::setdiff, list(...))}

# Related to strings ------------------


#' Extract the last n character(s) from a string
#' @description from https://stackoverflow.com/a/7963963/2116422
#' @param x a string
#' @param n number of characters to extract
#'
#' @return a vector of characters
#' @export
#'
substrRight <- function(x, n){ substr(x, nchar(x)-n+1, nchar(x))}

#' Remove empty spaces in strings
#'
#' @param x vector of characters
#'
#' @return the vector stripped of any space character
#' @export
#'
trimSpace<-function(x){
  gsub(" ","",as.character(x))
}

#' Capitalise the first character of a string
#'
#' @param x a vector of characters
#'
#' @return the vector of characters with the first one capitalised
#' @export
#'
simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

#' Capitalise the first character of a string
#'
#' @param x a vector of characters
#'
#' @return the vector of characters with the first one capitalised
#' @export
#'
simpleCap2<-function(x) {
  tryCatch(s<-strsplit(x, " ")[[1]],
           error = function(e) {s<-x})
}



#' Take a character string with parenthensis as a pattern and remove it from an another character string
#'
#' @param x character string that need to be modified
#' @param pattern character string to be removed
#'
#' @return a character string
#' @export
#'
cleanNames<-function(x,pattern){
  x<-gsub(paste(pattern,"("),'',x,fixed=TRUE)
  #x<-gsub(" (",'',x,fixed=TRUE)
  return(gsub("[)]",'',x))
}



# Other functions ------------------

#' For debugging purposes
#' Screen or report output
#' Look for a global variable DEBUG
#' @param x the object to be printed
#' @export
printDebug<-function(x){
  if(sharedEnv$debug)
    print(x)
}

#' Allow to change the value of debug
#'
#' @param x boolean
#'
#' @return changes the value of sharedEnv$debug
#' @export
#'
#' @examples setDebug(TRUE)
setDebug<-function(x){
  if(isTRUE(x)) assign("debug",x,sharedEnv)
  message(paste('debug set to',sharedEnv$debug))
}
assign("debug",FALSE,envir = sharedEnv)
