## from gplots 3.0.1.2  https://cran.r-project.org/web/packages/gplots/index.html
## original authors: Gregory R. Warnes, Ben Bolker, Lodewijk Bonebakker, Robert Gentleman, Wolfgang Huber Andy Liaw, Thomas Lumley, Martin Maechler, Arni Magnusson, Steffen Moeller, Marc Schwartz, Bill Venables
## added because package is orphaned at the moment
## License: GPL-2 (see https://www.gnu.org/licenses/old-licenses/gpl-2.0.html)
## $Id: heatmap.2.R 2103 2016-03-25 17:11:26Z warnes $

here <- function() {}

plot.dendrogram <- stats:::plot.dendrogram
environment(plot.dendrogram) <- environment(here)

plotNodeLimit <- stats:::plotNodeLimit
environment(plotNodeLimit) <- environment(here)

.memberDend <- stats:::.memberDend
environment(.memberDend) <- environment(here)

.midDend <- stats:::.midDend
environment(.midDend) <- environment(here)

unByteCode <- function(fun)
{
  FUN <- eval(parse(text=deparse(fun)))
  environment(FUN) <- environment(fun)
  FUN
}

plotNode <- unByteCode(stats:::plotNode)
environment(plotNode) <- environment(here)


#' Enhanced Heat Map (from gplots)
#' @description A heat map is a false color image (basically image(t(x))) with a dendrogram added to the left side and/or to the top. Typically, reordering of the rows and columns according to some set of values (row or column means) within the restrictions imposed by the dendrogram is carried out.
#' This heatmap provides a number of extensions to the standard R heatmap function.
#'
#' @param x numeric matrix of the values to be plotted.
#' @param Rowv determines if and how the row dendrogram should be reordered.	By default, it is TRUE, which implies dendrogram is computed and reordered based on row means. If NULL or FALSE, then no dendrogram is computed and no reordering is done. If a dendrogram, then it is used "as-is", ie without any reordering. If a vector of integers, then dendrogram is computed and reordered based on the order of the vector.
#' @param Colv determines if and how the column dendrogram should be reordered.	Has the options as the Rowv argument above and additionally when x is a square matrix, Colv="Rowv" means that columns should be treated identically to the rows.
#' @param distfun function used to compute the distance (dissimilarity) between both rows and columns. Defaults to dist.
#' @param hclustfun function used to compute the hierarchical clustering when Rowv or Colv are not dendrograms. Defaults to hclust.
#' @param dendrogram character string indicating whether to draw 'none', 'row', 'column' or 'both' dendrograms. Defaults to 'both'. However, if Rowv (or Colv) is FALSE or NULL and dendrogram is 'both', then a warning is issued and Rowv (or Colv) arguments are honoured.
#' @param reorderfun function(d, w) of dendrogram and weights for reordering the row and column dendrograms. The default uses stats{reorder.dendrogram}
#' @param symm logical indicating if x should be treated symmetrically; can only be true when x is a square matrix.
#' @param scale character indicating if the values should be centered and scaled in either the row direction or the column direction, or none. The default is "none".
#' @param na.rm logical indicating whether NA's should be removed.
#' @param revC logical indicating if the column order should be reversed for plotting, such that e.g., for the symmetric case, the symmetry axis is as usual.
#' @param add.expr expression that will be evaluated after the call to image. Can be used to add components to the plot.
#' @param breaks (optional) Either a numeric vector indicating the splitting points for binning x into colors, or a integer number of break points to be used, in which case the break points will be spaced equally between min(x) and max(x).
#' @param symbreaks Boolean indicating whether breaks should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise.
#' @param col colors used for the image. Defaults to heat colors (heat.colors).
#' @param colsep (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
#' @param rowsep (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
#' @param sepcolor (optional) vector of integers indicating which columns or rows should be separated from the preceding columns or rows by a narrow space of color sepcolor.
#' @param sepwidth (optional) Vector of length 2 giving the width (colsep) or height (rowsep) the separator box drawn by colsep and rowsep as a function of the width (colsep) or height (rowsep) of a cell. Defaults to c(0.05, 0.05)
#' @param cellnote (optional) matrix of character strings which will be placed within each color cell, e.g. p-value symbols.
#' @param notecex (optional) numeric scaling factor for cellnote items.
#' @param notecol (optional) character string specifying the color for cellnote text. Defaults to "cyan". Can be a matrix to personalise each cell.
#' @param na.color Color to use for missing value (NA). Defaults to the plot background color.
#' @param trace character string indicating whether a solid "trace" line should be drawn across 'row's or down 'column's, 'both' or 'none'. The distance of the line from the center of each color-cell is proportional to the size of the measurement. Defaults to 'column'.
#' @param tracecol character string giving the color for "trace" line. Defaults to "cyan".
#' @param hline Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if trace is 'row' or 'both'. Vertical lines are only drawn if trace 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
#' @param vline Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if trace is 'row' or 'both'. Vertical lines are only drawn if trace 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
#' @param linecol Vector of values within cells where a horizontal or vertical dotted line should be drawn. The color of the line is controlled by linecol. Horizontal lines are only plotted if trace is 'row' or 'both'. Vertical lines are only drawn if trace 'column' or 'both'. hline and vline default to the median of the breaks, linecol defaults to the value of tracecol.
#' @param margins numeric vector of length 2 containing the margins (see par(mar= *)) for column and row names, respectively.
#' @param ColSideColors (optional) character vector of length ncol(x) containing the color names for a horizontal side bar that may be used to annotate the columns of x.
#' @param RowSideColors (optional) character vector of length nrow(x) containing the color names for a vertical side bar that may be used to annotate the rows of x.
#' @param cexRow positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
#' @param cexCol positive numbers, used as cex.axis in for the row or column axis labeling. The defaults currently only use number of rows or columns, respectively.
#' @param labRow character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
#' @param labCol character vectors with row and column labels to use; these default to rownames(x) or colnames(x), respectively.
#' @param srtRow angle of row/column labels, in degrees from horizontal
#' @param srtCol angle of row/column labels, in degrees from horizontal
#' @param adjRow 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation).
#' @param adjCol 2-element vector giving the (left-right, top-bottom) justification of row/column labels (relative to the text orientation).
#' @param offsetRow Number of character-width spaces to place between row/column labels and the edge of the plotting region.
#' @param offsetCol Number of character-width spaces to place between row/column labels and the edge of the plotting region.
#' @param colRow color of row/column labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
#' @param colCol color of row/column labels, either a scalar to set the color of all labels the same, or a vector providing the colors of each label item
#' @param key logical indicating whether a color-key should be shown.
#' @param keysize numeric value indicating the size of the key
#' @param density.info character string indicating whether to superimpose a 'histogram', a 'density' plot, or no plot ('none') on the color-key.
#' @param denscol character string giving the color for the density display specified by density.info, defaults to the same value as tracecol.
#' @param symkey Boolean indicating whether the color key should be made symmetric about 0. Defaults to TRUE if the data includes negative values, and to FALSE otherwise.
#' @param densadj Numeric scaling value for tuning the kernel width when a density plot is drawn on the color key. (See the adjust parameter for the density function for details.) Defaults to 0.25.
#' @param key.title main title of the color key. If set to NA no title will be plotted.
#' @param key.xlab x axis label of the color key. If set to NA no label will be plotted.
#' @param key.ylab y axis label of the color key. If set to NA no label will be plotted.
#' @param key.xtickfun function computing tick location and labels for the xaxis of the color key. Returns a named list containing parameters that can be passed to axis. See examples.
#' @param key.ytickfun function computing tick location and labels for the y axis of the color key. Returns a named list containing parameters that can be passed to axis. See examples.
#' @param key.par graphical parameters for the color key. Named list that can be passed to par.
#' @param main main, x- and y-axis titles; defaults to none.
#' @param xlab main, x- and y-axis titles; defaults to none.
#' @param ylab main, x- and y-axis titles; defaults to none.
#' @param lmat visual layout: position matrix, column height, column width. See below for details
#' @param lhei visual layout: position matrix, column height, column width. See below for details
#' @param lwid visual layout: position matrix, column height, column width. See below for details
#' @param extrafun A function to be called after all other work. See examples.
#' @param ... additional arguments passed on to image
#'
#' @return invisibly, a list with components
#' rowInd	row index permutation vector as returned by order.dendrogram.
#' colInd	column index permutation vector.
#' call the matched call
#' rowMeans, rowSDs
#' mean and standard deviation of each row: only present if scale="row"
#' colMeans, colSDs
#' mean and standard deviation of each column: only present if scale="column"
#' carpet	reordered and scaled 'x' values used generate the main 'carpet'
#' rowDendrogram	row dendrogram, if present
#' colDendrogram	column dendrogram, if present
#' breaks	values used for color break points
#' col colors used
#' vline center-line value used for column trace, present only if trace="both" or trace="column"
#' hline	center-line value used for row trace, present only if trace="both" or trace="row"
#' colorTable	 A three-column data frame providing the lower and upper bound and color for each bin
#' layout	A named list containing the values used for lmat, lhei, and lwid.
#' @export
#'

heatmap.2 <- function (x,

                       ## dendrogram control
                       Rowv = TRUE,
                       Colv=if(symm)"Rowv" else TRUE,
                       distfun = dist,
                       hclustfun = hclust,
                       dendrogram = c("both","row","column","none"),
                       reorderfun = function(d, w) reorder(d, w),
                       symm = FALSE,

                       ## data scaling
                       scale = c("none","row", "column"),
                       na.rm=TRUE,

                       ## image plot
                       revC = identical(Colv, "Rowv"),
                       add.expr,

                       ## mapping data to colors
                       breaks,
                       symbreaks=any(x < 0, na.rm=TRUE) || scale!="none",

                       ## colors
                       col="heat.colors",

                       ## block sepration
                       colsep,
                       rowsep,
                       sepcolor="white",
                       sepwidth=c(0.05,0.05),

                       ## cell labeling
                       cellnote,
                       notecex=1.0,
                       notecol="cyan",
                       na.color=par("bg"),

                       ## level trace
                       trace=c("column","row","both","none"),
                       tracecol="cyan",
                       hline=median(breaks),
                       vline=median(breaks),
                       linecol=tracecol,

                       ## Row/Column Labeling
                       margins = c(5, 5),
                       ColSideColors,
                       RowSideColors,
                       cexRow = 0.2 + 1/log10(nr),
                       cexCol = 0.2 + 1/log10(nc),
                       labRow = NULL,
                       labCol = NULL,
                       srtRow = NULL,
                       srtCol = NULL,
                       adjRow = c(0,NA),
                       adjCol = c(NA,0),
                       offsetRow = 0.5,
                       offsetCol = 0.5,
                       colRow = NULL,
                       colCol = NULL,

                       ## color key + density info
                       key = TRUE,
                       keysize = 1.5,
                       density.info=c("histogram","density","none"),
                       denscol=tracecol,
                       symkey = any(x < 0, na.rm=TRUE) || symbreaks,
                       densadj = 0.25,
                       key.title = NULL,
                       key.xlab = NULL,
                       key.ylab = NULL,
                       key.xtickfun = NULL,
                       key.ytickfun = NULL,
                       key.par=list(),

                       ## plot labels
                       main = NULL,
                       xlab = NULL,
                       ylab = NULL,

                       ## plot layout
                       lmat = NULL,
                       lhei = NULL,
                       lwid = NULL,

                       ## extras
                       extrafun=NULL,
                       ...
)
{
  scale01 <- function(x, low=min(x), high=max(x) )
  {
    x <- (x-low)/(high - low)
    x
  }

  retval <- list()

  scale <- if(symm && missing(scale)) "none" else match.arg(scale)
  dendrogram <- match.arg(dendrogram)
  trace <- match.arg(trace)
  density.info <- match.arg(density.info)

  if(length(col)==1 && is.character(col) )
    col <- get(col, mode="function")

  if(!missing(breaks) && any(duplicated(breaks)) )
    stop("breaks may not contian duplicate values")

  if(!missing(breaks) && (scale!="none"))
    warning("Using scale=\"row\" or scale=\"column\" when breaks are",
            "specified can produce unpredictable results.",
            "Please consider using only one or the other.")

  if ( is.null(Rowv) || any(is.na(Rowv)) )
    Rowv <- FALSE
  if ( is.null(Colv) || any(is.na(Colv)) )
    Colv <- FALSE
  else if( all(Colv=="Rowv") )
    Colv <- Rowv


  if(length(di <- dim(x)) != 2 || !is.numeric(x))
    stop("`x' must be a numeric matrix")

  nr <- di[1]
  nc <- di[2]

  if(nr <= 1 || nc <= 1)
    stop("`x' must have at least 2 rows and 2 columns")

  if(!is.numeric(margins) || length(margins) != 2)
    stop("`margins' must be a numeric vector of length 2")

  if(missing(cellnote))
    cellnote <- matrix("", ncol=ncol(x), nrow=nrow(x))

  if(!inherits(Rowv, "dendrogram")) {
    ## Check if Rowv and dendrogram arguments are consistent
    if (
      (
        ( is.logical(Rowv) && !isTRUE(Rowv) )
        ||
        ( is.null(Rowv) )
      )
      &&
      ( dendrogram %in% c("both","row") )
    )
    {
      warning("Discrepancy: Rowv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting row dendogram.")

      if (dendrogram=="both")
        dendrogram <- "column"
      else
        dendrogram <- "none"

    }
  }

  if(!inherits(Colv, "dendrogram")) {
    ## Check if Colv and dendrogram arguments are consistent
    if (
      (
        (is.logical(Colv) && !isTRUE(Colv) )
        ||
        (is.null(Colv))
      )
      &&
      ( dendrogram %in% c("both","column")) )
    {
      warning("Discrepancy: Colv is FALSE, while dendrogram is `",
              dendrogram, "'. Omitting column dendogram.")

      if (dendrogram=="both")
        dendrogram <- "row"
      else
        dendrogram <- "none"

    }
  }


  ## by default order by row/col mean
  ## if(is.null(Rowv)) Rowv <- rowMeans(x, na.rm = na.rm)
  ## if(is.null(Colv)) Colv <- colMeans(x, na.rm = na.rm)

  ## get the dendrograms and reordering indices

  ## if( dendrogram %in% c("both","row") )
  ##  { ## dendrogram option is used *only* for display purposes
  if(inherits(Rowv, "dendrogram"))
  {
    ddr <- Rowv ## use Rowv 'as-is', when it is dendrogram
    rowInd <- order.dendrogram(ddr)
    if(length(rowInd)>nr || any(rowInd<1 | rowInd > nr ))
      stop("Rowv dendrogram doesn't match size of x")
    if (length(rowInd) < nr)
      nr <- length(rowInd)
  }
  else if (is.integer(Rowv))
  {
    ## Compute dendrogram and do reordering based on given vector
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)

    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Rowv))
  { ## If TRUE, compute dendrogram and do reordering based on rowMeans
    Rowv <- rowMeans(x, na.rm = na.rm)
    distr <- distfun(x)
    hcr <- hclustfun(distr)
    ddr <- as.dendrogram(hcr)
    ddr <- reorderfun(ddr, Rowv)

    rowInd <- order.dendrogram(ddr)
    if(nr != length(rowInd))
      stop("row dendrogram ordering gave index of wrong length")
  }
  else if(!isTRUE(Rowv))
  {
    rowInd <- nr:1
    ddr <- as.dendrogram(hclust(dist(diag(nr))))
  }
  else
  {
    rowInd <- nr:1
    ddr <- as.dendrogram(Rowv)
  }

  if(inherits(Colv, "dendrogram"))
  {
    ddc <- Colv ## use Colv 'as-is', when it is dendrogram
    colInd <- order.dendrogram(ddc)
    if(length(colInd)>nc || any(colInd<1 | colInd > nc ))
      stop("Colv dendrogram doesn't match size of x")
    if (length(colInd) < nc)
      nc <- length(colInd)
  }
  else if(identical(Colv, "Rowv")) {
    if(nr != nc)
      stop('Colv = "Rowv" but nrow(x) != ncol(x)')
    if(exists("ddr"))
    {
      ddc <- ddr
      colInd <- order.dendrogram(ddc)
    } else
      colInd <- rowInd
  } else if(is.integer(Colv))
  {## Compute dendrogram and do reordering based on given vector
    distc <- distfun(if(symm)x else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)

    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if (isTRUE(Colv))
  {## If TRUE, compute dendrogram and do reordering based on rowMeans
    Colv <- colMeans(x, na.rm = na.rm)
    distc <- distfun(if(symm)x else t(x))
    hcc <- hclustfun(distc)
    ddc <- as.dendrogram(hcc)
    ddc <- reorderfun(ddc, Colv)

    colInd <- order.dendrogram(ddc)
    if(nc != length(colInd))
      stop("column dendrogram ordering gave index of wrong length")
  }
  else if(!isTRUE(Colv))
  {
    colInd <- 1:nc
    ddc <- as.dendrogram(hclust(dist(diag(nc))))
  }
  else
  {
    colInd <- 1:nc
    ddc <- as.dendrogram(Colv)
  }

  retval$rowInd <- rowInd
  retval$colInd <- colInd
  retval$call <- match.call()


  ## reorder x & cellnote
  x <- x[rowInd, colInd]
  x.unscaled <- x
  cellnote <- cellnote[rowInd, colInd]

  if(is.null(labRow))
    labRow <- if(is.null(rownames(x))) (1:nr)[rowInd] else rownames(x)
  else
    labRow <- labRow[rowInd]

  if(is.null(labCol))
    labCol <- if(is.null(colnames(x))) (1:nc)[colInd] else colnames(x)
  else
    labCol <- labCol[colInd]

  if(!is.null(colRow))
    colRow <- colRow[rowInd]

  if(!is.null(colCol))
    colCol <- colCol[colInd]

  if(scale == "row") {
    retval$rowMeans <- rm <- rowMeans(x, na.rm = na.rm)
    x <- sweep(x, 1, rm)
    retval$rowSDs <-  sx <- apply(x, 1, sd, na.rm = na.rm)
    x <- sweep(x, 1, sx, "/")
  }
  else if(scale == "column") {
    retval$colMeans <- rm <- colMeans(x, na.rm = na.rm)
    x <- sweep(x, 2, rm)
    retval$colSDs <-  sx <- apply(x, 2, sd, na.rm = na.rm)
    x <- sweep(x, 2, sx, "/")
  }

  ## Set up breaks and force values outside the range into the endmost bins
  if(missing(breaks) || is.null(breaks) || length(breaks)<1 )
  {
    if( missing(col) ||  is.function(col) )
      breaks <- 16
    else
      breaks <- length(col)+1
  }

  if(length(breaks)==1)
  {
    if(!symbreaks)
      breaks <- seq( min(x, na.rm=na.rm), max(x,na.rm=na.rm), length=breaks)
    else
    {
      extreme <- max(abs(x), na.rm=TRUE)
      breaks <- seq( -extreme, extreme, length=breaks )
    }
  }

  nbr <- length(breaks)
  ncol <- length(breaks)-1

  if(class(col)=="function")
    col <- col(ncol)

  min.breaks <- min(breaks)
  max.breaks <- max(breaks)

  x[x<min.breaks] <- min.breaks
  x[x>max.breaks] <- max.breaks


  ## Calculate the plot layout
  if( missing(lhei) || is.null(lhei) )
    lhei <- c(keysize, 4)

  if( missing(lwid) || is.null(lwid) )
    lwid <- c(keysize, 4)

  if( missing(lmat) || is.null(lmat) )
  {
    lmat <- rbind(4:3, 2:1)

    if(!missing(ColSideColors)) { ## add middle row to layout
      if(!is.character(ColSideColors) || length(ColSideColors) != nc)
        stop("'ColSideColors' must be a character vector of length ncol(x)")
      lmat <- rbind(lmat[1,]+1, c(NA,1), lmat[2,]+1)
      lhei <- c(lhei[1], 0.2, lhei[2])
    }

    if(!missing(RowSideColors)) { ## add middle column to layout
      if(!is.character(RowSideColors) || length(RowSideColors) != nr)
        stop("'RowSideColors' must be a character vector of length nrow(x)")
      lmat <- cbind(lmat[,1]+1, c(rep(NA, nrow(lmat)-1), 1), lmat[,2]+1)
      lwid <- c(lwid[1], 0.2, lwid[2])
    }

    lmat[is.na(lmat)] <- 0
  }

  if(length(lhei) != nrow(lmat))
    stop("lhei must have length = nrow(lmat) = ", nrow(lmat))

  if(length(lwid) != ncol(lmat))
    stop("lwid must have length = ncol(lmat) =", ncol(lmat))

  ## Graphics `output' -----------------------

  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  layout(lmat, widths = lwid, heights = lhei, respect = FALSE)

  plot.index <- 1

  ## draw the side bars
  if(!missing(RowSideColors)) {
    par(mar = c(margins[1],0, 0,0.5))
    image(rbind(1:nr), col = RowSideColors[rowInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  if(!missing(ColSideColors)) {
    par(mar = c(0.5,0, 0,margins[2]))
    image(cbind(1:nc), col = ColSideColors[colInd], axes = FALSE)
    plot.index <- plot.index + 1
  }
  ## draw the main carpet
  par(mar = c(margins[1], 0, 0, margins[2]))
  #if(scale != "none" || !symm)
  #  {
  x <- t(x)
  cellnote <- t(cellnote)
  #  }
  if(revC)
  { ## x columns reversed
    iy <- nr:1
    if(exists("ddr"))
      ddr <- rev(ddr)
    x <- x[,iy]
    cellnote <- cellnote[,iy]
  }
  else iy <- 1:nr

  ## display the main carpet
  image(1:nc, 1:nr, x, xlim = 0.5+ c(0, nc), ylim = 0.5+ c(0, nr),
        axes = FALSE, xlab = "", ylab = "", col=col, breaks=breaks,
        ...)
  retval$carpet <- x
  if(exists("ddr"))
    retval$rowDendrogram <- ddr
  if(exists("ddc"))
    retval$colDendrogram <- ddc
  retval$breaks <- breaks
  retval$col <- col

  ## fill 'na' positions with na.color
  if(!invalid(na.color) & any(is.na(x)))
  {
    mmat <- ifelse(is.na(x), 1, NA)
    image(1:nc, 1:nr, mmat, axes = FALSE, xlab = "", ylab = "",
          col=na.color, add=TRUE)
  }

  ## add column labels
  if(is.null(srtCol) && is.null(colCol))
    axis(1,
         1:nc,
         labels= labCol,
         las= 2,
         line= -0.5 + offsetCol,
         tick= 0,
         cex.axis= cexCol,
         hadj=adjCol[1],
         padj=adjCol[2]
    )
  else
  {
    if(is.null(srtCol) || is.numeric(srtCol))
    {
      if(missing(adjCol) || is.null(adjCol))
        adjCol=c(1,NA)

      if(is.null(srtCol))
        srtCol <- 90

      xpd.orig <- par("xpd")
      par(xpd=NA)
      xpos <- axis(1, 1:nc, labels=rep("", nc), las=2, tick=0)
      text(x=xpos,
           y=par("usr")[3] - (1.0 + offsetCol) * strheight("M"),
           labels=labCol,
           ##pos=1,
           adj=adjCol,
           cex=cexCol,
           srt=srtCol,
           col=colCol
      )
      par(xpd=xpd.orig)
    }
    else
      warning("Invalid value for srtCol ignored.")
  }

  ## add row labels
  if(is.null(srtRow) && is.null(colRow))
  {
    axis(4,
         iy,
         labels=labRow,
         las=2,
         line=-0.5+offsetRow,
         tick=0,
         cex.axis=cexRow,
         hadj=adjRow[1],
         padj=adjRow[2]
    )
  }
  else
  {
    if(is.null(srtRow) || is.numeric(srtRow))
    {
      xpd.orig <- par("xpd")
      par(xpd=NA)
      ypos <- axis(4, iy, labels=rep("", nr), las=2, line= -0.5, tick=0)
      text(x=par("usr")[2] + (1.0 + offsetRow) * strwidth("M"),
           y=ypos,
           labels=labRow,
           adj=adjRow,
           cex=cexRow,
           srt=srtRow,
           col=colRow
      )
      par(xpd=xpd.orig)
    }
    else
      warning("Invalid value for srtRow ignored.")
  }



  ## add row and column headings (xlab, ylab)
  if(!is.null(xlab)) mtext(xlab, side = 1, line = margins[1] - 1.25)
  if(!is.null(ylab)) mtext(ylab, side = 4, line = margins[2] - 1.25)

  ## perform user-specified function
  if (!missing(add.expr))
    eval(substitute(add.expr))

  ## add 'background' colored spaces to visually separate sections
  if(!missing(colsep))
    for(csep in colsep)
      rect(xleft =csep+0.5,               ybottom=0,
           xright=csep+0.5+sepwidth[1],   ytop=ncol(x)+1,
           lty=1, lwd=1, col=sepcolor, border=sepcolor)

  if(!missing(rowsep))
    for(rsep in rowsep)
      rect(xleft =0,          ybottom= (ncol(x)+1-rsep)-0.5,
           xright=nrow(x)+1,  ytop   = (ncol(x)+1-rsep)-0.5 - sepwidth[2],
           lty=1, lwd=1, col=sepcolor, border=sepcolor)


  ## show traces
  min.scale <- min(breaks)
  max.scale <- max(breaks)
  x.scaled  <- scale01(t(x), min.scale, max.scale)

  if(trace %in% c("both","column") )
  {
    retval$vline <- vline
    vline.vals <- scale01(vline, min.scale, max.scale)
    for( i in 1:length(colInd) )
    {
      if(!is.null(vline))
      {
        abline(v=i-0.5 + vline.vals, col=linecol, lty=2)
      }
      xv <- rep(i, nrow(x.scaled)) + x.scaled[,i] - 0.5
      xv <- c(xv[1], xv)
      yv <- 1:length(xv)-0.5
      lines(x=xv, y=yv, lwd=1, col=tracecol, type="s")
    }
  }


  if(trace %in% c("both","row") )
  {
    retval$hline <- hline
    hline.vals <- scale01(hline, min.scale, max.scale)
    for( i in 1:length(rowInd) )
    {
      if(!is.null(hline))
      {
        abline(h=i - 0.5 + hline.vals, col=linecol, lty=2)
      }
      yv <- rep(i, ncol(x.scaled)) + x.scaled[i,] - 0.5
      yv <- rev(c(yv[1], yv))
      xv <- length(yv):1-0.5
      lines(x=xv, y=yv, lwd=1, col=tracecol, type="s")
    }
  }



  if(!missing(cellnote))
    text(x=c(row(cellnote)),
         y=c(col(cellnote)),
         labels=c(cellnote),
         col=c(notecol),
         cex=c(notecex))

  plot.index <- plot.index + 1

  ## increment plot.index and then do
  ##   latout_set( lmat, plot.index )
  ## to set to the correct plot region, instead of
  ## relying on plot.new().

  ## the two dendrograms :
  par(mar = c(margins[1], 0, 0, 0))
  if( dendrogram %in% c("both","row") )
  {
    flag <- try(
      plot.dendrogram(ddr, horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none")
    )
    if("try-error" %in% class(flag))
    {
      cond <- attr(flag, "condition")
      if(!is.null(cond) && conditionMessage(cond)=="evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop('Row dendrogram too deeply nested, recursion limit exceeded.  Try increasing option("expressions"=...).')
    }
  }
  else
    plot.new()

  par(mar = c(0, 0, if(!is.null(main)) 5 else 0, margins[2]))

  if( dendrogram %in% c("both","column") )
  {
    flag <- try(
      plot.dendrogram(ddc, axes = FALSE, xaxs = "i", leaflab = "none")
    )
    if("try-error" %in% class(flag))
    {
      cond <- attr(flag, "condition")
      if(!is.null(cond) && conditionMessage(cond)=="evaluation nested too deeply: infinite recursion / options(expressions=)?")
        stop('Column dendrogram too deeply nested, recursion limit exceeded.  Try increasing option("expressions"=...).')
    }
  }
  else
    plot.new()

  ## title
  if(!is.null(main)) title(main, cex.main = 1.5*op[["cex.main"]])

  ## Add the color-key
  if(key)
  {
    mar <- c(5, 4, 2, 1)
    if (!is.null(key.xlab) && is.na(key.xlab))
      mar[1] <- 2
    if (!is.null(key.ylab) && is.na(key.ylab))
      mar[2] <- 2
    if (!is.null(key.title) && is.na(key.title))
      mar[3] <- 1
    par(mar = mar, cex=0.75, mgp=c(2, 1, 0))
    if (length(key.par) > 0)
      do.call(par, key.par)
    tmpbreaks <- breaks

    if(symkey)
    {
      max.raw <- max(abs(c(x,breaks)),na.rm=TRUE)
      min.raw <- -max.raw
      tmpbreaks[1] <- -max(abs(x), na.rm=TRUE)
      tmpbreaks[length(tmpbreaks)] <- max(abs(x), na.rm=TRUE)
    }
    else
    {
      min.raw <- min.breaks
      max.raw <- max.breaks
    }

    z <- seq(min.raw, max.raw, by=min(diff(breaks)/100))
    image(z=matrix(z, ncol=1),
          col=col, breaks=tmpbreaks,
          xaxt="n", yaxt="n")

    par(usr=c(0,1,0,1))
    if (is.null(key.xtickfun)) {
      lv <- pretty(breaks)
      xv <- scale01(as.numeric(lv), min.raw, max.raw)
      xargs <- list(at=xv, labels=lv)
    } else {
      xargs <- key.xtickfun()
    }
    xargs$side <- 1
    do.call(axis, xargs)
    if (is.null(key.xlab)) {
      if(scale=="row")
        key.xlab <- "Row Z-Score"
      else if(scale=="column")
        key.xlab <- "Column Z-Score"
      else
        key.xlab <- "Value"
    }
    if (!is.na(key.xlab)) {
      mtext(side=1, key.xlab, line=par("mgp")[1], padj=0.5, cex=par("cex") * par("cex.lab"))
    }

    if(density.info=="density")
    {
      dens <- density(x, adjust=densadj, na.rm=TRUE,
                      from=min.scale, to=max.scale)
      omit <- dens$x < min(breaks) | dens$x > max(breaks)
      dens$x <- dens$x[!omit]
      dens$y <- dens$y[!omit]
      dens$x <- scale01(dens$x, min.raw, max.raw)
      lines(dens$x, dens$y / max(dens$y) * 0.95, col=denscol, lwd=1)
      if (is.null(key.ytickfun)) {
        yargs <- list(at=pretty(dens$y)/max(dens$y) * 0.95, labels=pretty(dens$y))
      } else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Density Plot"
      if (!is.na(key.title))
        title(key.title)
      par(cex=0.5)
      if (is.null(key.ylab))
        key.ylab <- "Density"
      if (!is.na(key.ylab))
        mtext(side=2,key.ylab, line=par("mgp")[1], padj=0.5, cex=par("cex") * par("cex.lab"))
    }
    else if(density.info=="histogram")
    {
      h <- hist(x, plot=FALSE, breaks=breaks)
      hx <- scale01(breaks, min.raw, max.raw)
      hy <- c(h$counts, h$counts[length(h$counts)])
      lines(hx, hy/max(hy)*0.95, lwd=1, type="s", col=denscol)
      if (is.null(key.ytickfun)) {
        yargs <- list(at=pretty(hy)/max(hy) * 0.95, labels=pretty(hy))
      } else {
        yargs <- key.ytickfun()
      }
      yargs$side <- 2
      do.call(axis, yargs)
      if (is.null(key.title))
        key.title <- "Color Key\nand Histogram"
      if (!is.na(key.title))
        title(key.title)
      par(cex=0.5)
      if (is.null(key.ylab))
        key.ylab <- "Count"
      if (!is.na(key.ylab))
        mtext(side=2,key.ylab, line=par("mgp")[1], padj=0.5, cex=par("cex") * par("cex.lab"))
    }
    else
      if (is.null(key.title))
        title("Color Key")

    if(trace %in% c("both","column") )
    {
      vline.vals <- scale01(vline, min.raw, max.raw)
      if(!is.null(vline))
      {
        abline(v=vline.vals, col=linecol, lty=2)
      }
    }


    if(trace %in% c("both","row") )
    {
      hline.vals <- scale01(hline, min.raw, max.raw)
      if(!is.null(hline))
      {
        abline(v=hline.vals, col=linecol, lty=2)

      }
    }

  }
  else
  {
    par(mar=c(0, 0, 0, 0))
    plot.new()
  }
  ## Create a table showing how colors match to (transformed) data ranges
  retval$colorTable <- data.frame(
    low=retval$breaks[-length(retval$breaks)],
    high=retval$breaks[-1],
    color=retval$col
  )

  # Store layout information, suggested by Jenny Drnevich
  retval$layout <- list(lmat = lmat,
                        lhei = lhei,
                        lwid = lwid
  )


  ## If user has provided an extra function, call it.
  if(!is.null(extrafun))
    extrafun()

  invisible( retval )
}

# $Id: colorpanel.R 736 2005-11-18 00:16:28Z warnes $

#' Generate a smoothly varying set of colors
#' @description colorpanel generate a set of colors that varies smoothly. redgreen, greenred, bluered, and redblue generate red-black-green, green-black-red, red-white-blue, and blue-white-red colorbars, respectively. colors
#'
#' @param n Desired number of color elements in the panel.
#' @param low Color to use for the lowest value.
#' @param mid Color to use for the  middle value. May be ommited.
#' @param high Color to use for the highest values.
#'
#' @return Vector of HTML-style RGB colors.
#' @export
#'

colorpanel <- function(n,low,mid,high)
{
  if(missing(mid) || missing(high) )
  {
    ## convert to rgb
    low <- col2rgb(low)
    if(missing(high))
      high <- col2rgb(mid)
    else
      high <- col2rgb(high)

    red    <- seq(low[1,1], high[1,1], length=n)/255
    green  <- seq(low[3,1], high[3,1], length=n)/255
    blue   <- seq(low[2,1], high[2,1], length=n)/255
  }
  else # use a center color
  {
    isodd <- odd(n)
    if(isodd)
    {
      n <- n+1
    }

    ## convert to rgb
    low <- col2rgb(low)
    mid <- col2rgb(mid)
    high <- col2rgb(high)

    ## determine length of each component
    lower <- floor(n/2)
    upper <- n - lower

    red  <- c(
      seq(low[1,1], mid [1,1], length=lower),
      seq(mid[1,1], high[1,1], length=upper)
    )/255

    green <- c(
      seq(low[3,1], mid [3,1], length=lower),
      seq(mid[3,1], high[3,1], length=upper)
    )/255

    blue <- c(
      seq(low[2,1], mid [2,1], length=lower),
      seq(mid[2,1], high[2,1], length=upper)
    )/255

    if(isodd)
    {
      red   <- red  [-(lower+1)]
      green <- green[-(lower+1)]
      blue  <- blue [-(lower+1)]
    }
  }

  rgb(red,blue,green)
}




#' Generate red-to-green colorscale
#'
#' @param n Desired number of color elements in the panel.
#'
#' @return Vector of HTML-style RGB colors.
#' @export
#'
redgreen <- function(n) {colorpanel(n, 'red', 'black', 'green')}

#' Generate green-to-red colorscale
#'
#' @param n Desired number of color elements in the panel.
#'
#' @return Vector of HTML-style RGB colors.
#' @export
#'
greenred <- function(n) {colorpanel(n, 'green', 'black', 'red' )}

#' Generate blue-to-red colorscale
#'
#' @param n  Desired number of color elements in the panel.
#'
#' @return Vector of HTML-style RGB colors.
#' @export
#'
bluered  <- function(n) {colorpanel(n, 'blue','white','red')}

#' Generate red-to-blue colorscale
#'
#' @param n Desired number of color elements in the panel.
#'
#' @return Vector of HTML-style RGB colors.
#' @export
#'
redblue  <- function(n) {colorpanel(n, 'red','white','blue')}





