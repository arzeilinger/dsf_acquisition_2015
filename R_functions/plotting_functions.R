#### Plotting functions from Daniel Turek

#### samplesPlot function
samplesPlot <- function(samples, var=colnames(samples), ind=NULL, burnin=NULL, width=7, height=4, legend=TRUE, legend.location='topright', traceplot=TRUE, densityplot=TRUE, file=NULL) {
  if(!is.null(file)) tiff(file, width=width, height=height, res = 300, compression = "lzw", units = "cm") else
    if(inherits(try(knitr::opts_chunk$get('dev'), silent=TRUE), 'try-error') || is.null(knitr::opts_chunk$get('dev')))   ## if called from Rmarkdown/knitr
      dev.new(height=height, width=width)
  par.save <- par(no.readonly = TRUE)
  par(mfrow=c(1,traceplot+densityplot), cex=0.7, cex.main=1.5, cex.axis=0.9, cex.lab=1.5, lab=c(3,3,7), mgp=c(1.5,0.4,0), mar=c(3,3,2,0.6), oma=c(0,0,0,0), tcl=-0.3, bty='l')
  ## process samples
  var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', var))   ## add \\ before any '[' or ']' appearing in var
  var <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
  samples <- samples[, var, drop=FALSE]
  if(!is.null(ind) && !is.null(burnin)) stop('only specify either ind or burnin')
  if(!is.null(ind))     samples <- samples[ind, , drop=FALSE]
  if(!is.null(burnin))  samples <- samples[(burnin+1):dim(samples)[1], , drop=FALSE]
  nparam <- ncol(samples)
  rng <- range(samples)
  if(!traceplot & !densityplot) stop('both traceplot and densityplot are false')
  if(traceplot) {  ## traceplot
    plot(1:nrow(samples), ylim=rng, type='n', main='Traceplots', xlab='', ylab='')
    for(i in 1:nparam)
      lines(samples[,i], col=rainbow(nparam, alpha=0.75)[i])
    if(legend & !densityplot & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
      legend(legend=dimnames(samples)[[2]], fill=rainbow(nparam, alpha=0.5), bty='n', x=legend.location)
  }  ## finish traceplot
  if(densityplot) {  ## denstyplot
    xMin <- xMax <- yMax <- NULL
    for(i in 1:nparam) {
      d <- density(samples[,i])
      xMin <- min(xMin,d$x); xMax <- max(xMax,d$x); yMax <- max(yMax, d$y) }
    plot(1, xlim=c(xMin,xMax), ylim=c(0,yMax), type='n', main='', xlab = "Posterior probability", ylab = "Density", yaxt='n')
    for(i in 1:nparam)
      polygon(density(samples[,i]), col=grey.colors(nparam, alpha=0.5)[i], border="black",
              main = "")
      # To keep transparent border
      # polygon(density(samples[,i]), col=grey.colors(nparam, alpha=0.5)[i], border=grey.colors(nparam, alpha=0.5)[i])
    if(legend & !is.null(dimnames(samples)) & is.character(dimnames(samples)[[2]]))
      legend(legend=dimnames(samples)[[2]], fill=grey.colors(nparam, alpha=0.5), bty='n', x=legend.location)
  }  ## finish densityplot
  if(!is.null(file)) dev.off()
  invisible(par(par.save))
}


#### chainsPlot function
chainsPlot <- function(samplesList, var=NULL, nrows=3, width=7, height=min(1+3*nrows,7), legend=!is.null(names(samplesList)), legend.location='topright', jitter=1, buffer.right=0, buffer.left=0, cex=1, file=NULL) {
  if(!is.null(file)) pdf(file, width=width, height=height) else
    if(inherits(try(knitr::opts_chunk$get('dev'), silent=TRUE), 'try-error') || is.null(knitr::opts_chunk$get('dev')))   ## if called from Rmarkdown/knitr
      dev.new(height=height, width=width)
  par.save <- par(no.readonly = TRUE)
  par(mfrow=c(nrows,1), oma=c(3,1,1,1), mar=c(4,1,0,1), mgp=c(3,0.5,0))
  if(!(class(samplesList) %in% c('list', 'mcmc.list'))) samplesList <- list(samplesList)
  if(!is.null(var)) samplesList <- lapply(samplesList, function(samples) {
    var <- gsub('\\[', '\\\\\\[', gsub('\\]', '\\\\\\]', var))   ## add \\ before any '[' or ']' appearing in var
    theseVar <- unlist(lapply(var, function(n) grep(paste0('^', n,'(\\[.+\\])?$'), colnames(samples), value=TRUE)))  ## expanded any indexing
    samples[, theseVar, drop=FALSE]
  })
  chainParamNamesList <- lapply(samplesList, function(s) colnames(s))
  nChains <- length(samplesList)
  paramNamesAll <- unique(unlist(lapply(samplesList, function(s) colnames(s))))
  nParamsAll <- length(paramNamesAll)
  cols <- rainbow(nChains)
  ## construct 3D summary array:
  summary <- array(as.numeric(NA), dim = c(nChains, 3, nParamsAll))
  if(!is.null(names(samplesList))) dimnames(summary)[[1]] <- names(samplesList)
  dimnames(summary)[[2]] <- c('mean','low','upp')
  dimnames(summary)[[3]] <- paramNamesAll
  for(iChain in 1:nChains) {
    theseSamples <- samplesList[[iChain]]
    thisSummary <- rbind(mean = apply(theseSamples, 2, mean),
                         low  = apply(theseSamples, 2, function(x) quantile(x, 0.025)),
                         upp  = apply(theseSamples, 2, function(x) quantile(x, 0.975)))
    summary[iChain,c('mean','low','upp'),colnames(thisSummary)] <- thisSummary
  }
  nParamsPerRow <- ceiling(nParamsAll/nrows)
  sq <- if(nChains==1) 0 else seq(-1,1,length=nChains)
  scale <- width/nParamsPerRow * jitter * 0.1  ## adjust jitter scale factor at end
  for(iRow in 1:nrows) {
    rowParamInd <- (1+(iRow-1)*nParamsPerRow) : ifelse(iRow==nrows,nParamsAll,iRow*nParamsPerRow)
    nRowParams <- length(rowParamInd)
    rowParamNames <- paramNamesAll[rowParamInd]
    xs <- 1:nRowParams
    names(xs) <- rowParamNames
    ylim <- range(summary[,c('low','upp'),rowParamNames], na.rm=TRUE)
    plot(x=-100, y=0, xlim=c(1-buffer.left,nParamsPerRow+buffer.right), ylim=ylim, xaxt='n', ylab='', xlab='', tcl=-0.3, cex.axis=cex)
    axis(1, at=1:nRowParams, labels=FALSE, tcl=-0.3)
    text(x=1:nRowParams, y=par()$usr[3]-0.1*(par()$usr[4]-par()$usr[3]), labels=rowParamNames, srt=45, adj=1, xpd=TRUE, cex=0.9*cex)
    for(iChain in 1:nChains) {
      ps <- intersect(rowParamNames, chainParamNamesList[[iChain]])
      xsJittered <- xs + sq[iChain]*scale
      points(x=xsJittered[ps], y=summary[iChain,'mean',ps], pch=16, col=cols[iChain])
      segments(x0=xsJittered[ps], y0=summary[iChain,'low',ps], y1=summary[iChain,'upp',ps], lwd=1, col=cols[iChain])
    }
    if(legend) legend(legend.location, legend=names(samplesList), pch=16, col=cols, cex=cex)
  }
  if(!is.null(file)) dev.off()
  invisible(par(par.save))
}

