require(flowCore)
require(sfsmisc)
require(changepoint)
require(bit)
require(grid)
Sys.setlocale('LC_ALL','C') ## Occasional multibyte string issue

geo.mean <- function(vec){
  return(exp(mean(log(vec))))
}

make_pops <- function(dF, cutoff, params, markers){
  dF <- dF[,params]
  cnames <- colnames(dF) 
  if (cutoff == "median"){ cutoff <- apply(dF, 2, function(x){ quantile(x, 0.5) })}
  else if (cutoff < 1){
    cutoff <- apply(dF, 2, function(x, v) { quantile(x, v) }, v=cutoff)
  }
  else { cutoff <- rep(cutoff, length(params)) }
  d2 <- do.call(cbind, lapply(1:ncol(dF), function(i, a, b){
    x <- a[,i]
    cutoff <- b[i]
    if ( (length(x[which(x > cutoff)])/length(x)) > 0.9){
      cutoff <- median(x)
    }
    x[which(x <= cutoff)] <- 0
    x[which(x != 0)] <- 1
    return(x)
  }, a=dF, b=cutoff))

  dfbit <- apply(d2, 1, as.bit)
  aphbit <- unique(dfbit)
  idx <- lapply(c(1:length(aphbit)), function(x, y, i){
    return(which(x == y[[i]]))
  }, x=dfbit, y=aphbit)

  lengths <- lapply(idx, length)
  return(list("aphbit"=aphbit, "lengths"=lengths, "idx"=idx))
}

name_maker <- function(aphbit, markers, full=TRUE){
  v <- lapply(aphbit, function(x, d){
    if (x == 0) { return(rep.int(0, times=d)) }
    else {
      dig <- t(digitsBase(x, base=2))[1:length(markers)]
      if (length(dig) == length(markers)){ return(dig) }
      else {
        dig <- c(rep.int(0, times=(length(markers) - length(dig))), dig)
        return(dig)
      }
     } },
    d=length(aphbit))

  newv <- lapply(v, function(w, x, full){
    ids <- !is.na(x)
    x[which(x == 1)] <- "+"
    x[which(x == 0)] <- "-"
    x <- lapply(seq_along(x), function(i, y, z){
      return(paste(z[i], y[i], sep=""))
    }, y=x, z=w)
    if (full == FALSE) { x <- x[ids] }
    x <- paste(x, collapse="")
    return(x)
  }, w=markers, full=full)
  return(newv)
}

get_pops <- function(dF, cutoff, params, bins, nCellCutoff, markers, nstable){
  poplist <- make_pops(dF, cutoff, params, markers)
  curr <- length(which(poplist$lengths >= nCellCutoff))
  ## try to find as many pops as possible
  if (length(which(poplist$lengths >= nCellCutoff)) < nstable){
    cut_grid <- rev(seq(0.05, 0.5, by=0.05))
    while (curr < 5 & length(cut_grid) >= 1){
      poplist <- make_pops(dF, cut_grid[1], params, markers)
      curr <- length(which(poplist$lengths >= nCellCutoff))
      cut_grid <- cut_grid[-1]
    }
    ## logic for seemingly redundant code:
    ## if the conditional evaluates and is true there are few pops
    ## if the conditional doesn't evaluate there are definitely few pops
    tryCatch({
        if(cut_grid == 0.05 & curr < nstable){
            warning(paste0("Less than ", nstable," populations identified by flowClean"))
        }
    }, error = function(e){
        warning(paste0("Less than ", nstable," populations identified by flowClean"))
    })
  }
  aphbit <- poplist$aphbit
  lengths <- poplist$lengths
  idx <- poplist$idx
  pops <- name_maker(aphbit, markers)

  counts <- lapply(1:length(idx), function(i, x, y){
    cc <- lapply(1:length(bins), function(j, w, v) {
      return(length(which(w %in% v[[j]])))
    }, w=x[[i]], v=bins)
    return(cc)
  }, x=idx, y=bins)

  popdf <- do.call(rbind, counts)
  perdf <- apply(popdf, 2, function(x){ N <- unlist(x); return(N/sum(N)) })
  perdf <- cbind(perdf, as.double(lengths))
  rownames(perdf) <- pops
  good.idx <- which(lengths >= nCellCutoff)
  perdf.trim <- perdf[good.idx,1:length(bins)]
  return(list("full"=perdf, "trim"=perdf.trim, "pops"=popdf))
}

clean <- function(fF, vectMarkers, filePrefixWithDir, ext, binSize=0.01, nCellCutoff=500,
                  announce=TRUE, cutoff="median", diagnostic=FALSE, fcMax=1.3,
                  returnVector=FALSE, nstable=5){

  if (dim(exprs(fF))[1] < 10000){
      warning("Too few cells in FCS for flowClean.")
      GoodVsBad <- rep.int(0, times=nrow(exprs(fF)))
      outFCS <- makeFCS(fF, GoodVsBad, filePrefixWithDir, 0, nCellCutoff, ext,
                        stablePops=NULL)
      return(outFCS)
  }

  markers <- parameters(fF)$name
  markers <- as.vector(markers)

  numbins <- ceiling(1/binSize)
  time.id <- grep("time", colnames(exprs(fF)), ignore.case = TRUE)
  if (length(time.id) > 0) {
      time <- exprs(fF)[, time.id]
      ## Have to deal with this here
      ## if (mean(time) == time[1]) {
      ##     time <- 1:nrow(exprs(fF))
      ## }
  }
  else {
      warning("No Time Parameter Detected")
      GoodVsBad <- rep.int(0, times=nrow(exprs(fF)))
      outFCS <- makeFCS(fF, GoodVsBad, filePrefixWithDir, 0, nCellCutoff, ext,
                        stablePops=NULL)
      return(outFCS)
  }
  # make sure time starts at 0
  if (min(time) > 0){ time <- time - min(time) }
  numOfEvents <- length(time)
  lTime <- time[numOfEvents]
  stepB <- lTime * binSize
  bins <- lapply(c(1:numbins), function(i, x, y, z){
    vec <- which(x >= ((i - 1) * y) & (x < i * y))
    if (i == z){ vec <- c(vec, which(x >= (i * y))) }
    return(vec)
  }, x=time, y=stepB, z=numbins )
  binVector <- unlist(lapply(c(1:numbins), function(i, x){
    rep(i, length(unlist(x[[i]]))) }, x=bins))
  ### NOTE: implicitly assuming vectMarkers is numeric
  out <- get_pops(exprs(fF), cutoff, params=vectMarkers, bins=bins,
                  nCellCutoff, markers[vectMarkers], nstable)
  full <- out$full
  pops <- out$pops
  out <- out$trim
  dxVector <- binVector
  ## missing cells really screw everything up
  pops <- apply(pops, c(1, 2), as.numeric)
  sums <- apply(pops, 2, sum)
  weirds <- lapply(c(37, 42, 51), function(xx){
      weird <- getUnlikely(sums, xx)
  })
  weird <- sort(unique(unlist(weirds)))
  if (length(weird) > 0){
      out <- out[,-weird]
  }
  out <- cen.log.ratio(out)
  norms <- lp(out)
  ## was previously penalty=AIC, but with recent updates this works
  ## better/does what AIC used to
  ## what is the indexing here?  
  pts <- cpt.mean(norms, method="PELT", penalty="Manual", pen.value=1)
  bad <- getBad(pts, fcMax)

  ## what kind of bad do we report?
  if (!is.null(bad) | length(weird) > 0){
    if (length(weird) > 0){
        bad <- unique(c(unlist(fix.weird(bad, weird, numbins)), weird))
    }
    else if (bad[length(bad)] != numbins){
        ## changepoint library frequently off by 1
        bad <- c(min(bad), unique(bad+1))
    }
    bad <- sort(bad)
    dxVector[which(dxVector %in% bad)] <- runif(length(which(dxVector %in% bad)),
                                                min=10000, max=20000)
    GoodVsBad <- as.numeric(dxVector)
    if (returnVector == TRUE){ return(GoodVsBad) }
    if (diagnostic){
      png(paste(filePrefixWithDir,sep=".", numbins, nCellCutoff,
                "clr_percent_plot", "png"), type="cairo",
          height=1000, width=1000)
      diagnosticPlot(out,"CLR", bad)
      dev.off()
    }

    outFCS <- makeFCS(fF, GoodVsBad, filePrefixWithDir, numbins, nCellCutoff,
                      ext, stablePops=out)
    if (announce){
      print(paste("flowClean has identified problems in ",
                  description(fF)$FILENAME, " with ", toString(bad),  ".", sep=""))
    }
    return(outFCS)
  }
  else{
    if (announce){
        print(paste("flowClean detected no problems in ",
                    description(fF)$FILENAME, ".", sep=""))
    }
    if (diagnostic){
      png(paste(filePrefixWithDir,sep=".", numbins, nCellCutoff,
                "clr_percent_plot", "png"), type="cairo",
          height=1000, width=1000)
      diagnosticPlot(out,"CLR", bad)
      dev.off()
   }

    GoodVsBad <- as.numeric(dxVector)
    if (returnVector == TRUE){ return(GoodVsBad) }
    outFCS <- makeFCS(fF, GoodVsBad, filePrefixWithDir, numbins, nCellCutoff,
                      ext, stablePops=out)
    return(outFCS)
  } 
}

makeSeries <- function(vec){
    out <- list()
    out[[1]] <-  holder <- vec[1]
    counter <- 1
    devnull <- lapply(diff(vec), function(xx){
        holder <<- holder + xx
        if (xx != 1){
            counter <<- counter + 1
            out[[counter]] <<- holder
        }
        else {
            out[[counter]] <<- c(out[[counter]], holder)
        }
    })
    out
}

fix.weird <- function(bad, weird, maxBin){
    bad <- c(min(bad)-1, bad)
    if(is.null(bad)){ return(weird) }
    ser.b <- makeSeries(bad)
    ser.w <- makeSeries(weird)
    ## fix persistent off-by-one error
    ## need to include minimum bad
    ser.b <- lapply(ser.b, function(xx){
        ## with gaps in the data, we also miss the start of problems
        if (min(xx) > 1){ xx <- c(min(xx) - 1, xx) }
        if (max(xx) < maxBin){ return(xx + 1) }
        else{ xx }
    })
    ## let's find the max per weird series and place it before
    ## a bad series if possible
    ## creates a vector whose indices are those of
    ## 'weird series' and whose entries index 'bad series'
    shifts <- sapply(sapply(ser.w, max), function(xx){
        which(sapply(ser.b, min) > xx)[1]
    })
    ## occasionally the weird bin is the last bin
    if (!is.na(shifts)){
      devnull <- sapply(1:length(shifts[!is.na(shifts)]), function(ii){
        bad.id <- shifts[ii]
        bad.offset <- length(ser.w[[ii]])
        ser.b[[bad.id]] <<- ser.b[[bad.id]] + bad.offset
      })
    }
    ## problem is not here!
    ser.b
}

loglik <- function(data, l){
    n <- length(data)
    sum(data, na.rm=TRUE) * log(l) - n*l
}

getUnlikely <- function(sums, seed){
   set.seed(seed)
   idx <- sample(1:length(sums))
   names(sums) <- paste0("V", 1:length(sums))
   sums <- sums[idx]
   cll <- sapply(1:length(sums), function(ii){
       if (ii == 1){return(loglik(sums[ii], sums[ii]))}
       l <- loglik(sums[1:(ii-1)], mean(sums[1:(ii-1)]))
       lnew <- loglik(sums[1:ii], mean(sums[1:ii]))
       return(lnew - l)
   })
   weird <- which(cll <= 0)
   weird <- as.numeric(gsub("V", "", names(sums)[weird]))
   weird
}


makeFCS <- function(fF, GoodVsBad, filePrefixWithDir, numbins, nCellCutoff, ext, stablePops){
  ex <- exprs(fF)
  rs <- attr(exprs(fF), "ranges")
  rs <- c(rs, rs[1])
  ex <- cbind(ex, GoodVsBad)
  attr(ex, "ranges") <- rs
  NN <- as.numeric(keyword(fF)["$PAR"]) + 1
  names(dimnames(ex)[[2]]) <- sprintf("$P%sN", 1:NN)
  pnr <- paste0("$P", NN, "R")
  pnb <- paste0("$P", NN, "B")
  pne <- paste0("$P", NN, "E")
  pnn <- paste0("$P", NN, "N")
  pns <- paste0("$P", NN, "S")
  flowCorePnRmax <- paste0("flowCore_$P", NN, "Rmax")
  flowCorePnRmin <- paste0("flowCore_$P", NN, "Rmin")

  o <- parameters(fF)@data
  o[length(o[,1]) + 1,] <- c("GoodVsBad", "GoodVsBad", as.numeric(keyword(fF)$`$P1R`), 0, as.numeric(keyword(fF)$`$P1R`) - 1)
  outFCS <- new("flowFrame", exprs=ex, parameters=new("AnnotatedDataFrame",o), keyword1=keyword(fF))
  keyword1(outFCS)$FILENAME <- paste(filePrefixWithDir,sep=".", numbins, nCellCutoff, "revised", ext)
  keyword1(outFCS)[pnr] <- max(20000, keyword1(outFCS)$`$P1R`)
  keyword1(outFCS)[pnb] <- keyword1(outFCS)$`$P1B`
  keyword1(outFCS)[pne] <- "0,0"
  keyword1(outFCS)[pnn] <- "GoodVsBad"
  keyword1(outFCS)[pns] <- "GoodVsBad"
  keyword1(outFCS)$`$PAR` <- NN
  keyword1(outFCS)$`StablePops` <- nrow(stablePops)
  keyword1(outFCS)$`nBins` <- numbins
  keyword1(outFCS)$`nCellCutoff` <- nCellCutoff
  keyword1(outFCS)[flowCorePnRmax] <- max(20000, keyword1(outFCS)$`flowCore_$P1Rmax`)
  keyword1(outFCS)[flowCorePnRmin] <- 0
  parameters(outFCS)@data$range <- as.numeric(parameters(outFCS)@data$range)
  parameters(outFCS)@data$minRange <- as.numeric(parameters(outFCS)@data$minRange)
  parameters(outFCS)@data$maxRange <- as.numeric(parameters(outFCS)@data$maxRange)
  outFCS
}  

diagnosticPlot <- function(dF, ylab, bad){
  bins <- ncol(dF)
  pops <- nrow(dF)
  xx <- 1:bins
  plot(xx, dF[1,],ylim=c(min(dF), max(dF)),xlab="Time Bins", ylab=ylab, col=ifelse(xx %in% bad, 2, 1), pch=20)
  for (i in 2:pops){
    points(1:bins,dF[i,], col=ifelse(xx %in% bad, 2, 1), pch=20)
  }
  abline(h=0, col=4)
}

cen.log.ratio <- function(dF, minim=1e-7){
  #element-wise statistic
  sums <- unlist(apply(dF, 1, sum))
  ta <- minim/max(sums)
  n <- nrow(dF)

  #element-wise statistic
  m <- max(unlist(apply(dF, 1, function(x){return(length(x[which(x == 0)]))})))
  if (n == m) { n <- n+1 }   
  d <- (n^2*ta)/((m+1)*(n-m))
  Ts <- (d*m*(m+1))/(n^2)

  #modified Aitchison of Fry et al
  rev_df <- apply(dF, c(1,2), function(x){
    if (x == 0){ x <- ta }
    else { x <- x - (x*Ts) }; return(x)
  })
  rev_df <- apply(rev_df, 2, function(x) {return(log(x/geo.mean(x)))})

  return(rev_df)
}

getCPTs <- function(cpt){
    methodFound <- slotFound <- FALSE
    out <- NULL
    if (length(cpts(cpt)) > 0) {
        methodFound <- TRUE
        out <- cpts(cpt)
    }
    else if (length(cpt@cpts)){
        slotFound <- TRUE
        out <- cpt@cpts
    }
    if (methodFound | slotFound){ return(out) }
    else { return(NULL) }
}

getBad <- function(cpt, k=1.3){
    pts <- getCPTs(cpt)
    if (length(pts) == 0) {
        bad <- NULL
        return(bad)
    }
    ps <- ps2 <- pts
    ms <- param.est(cpt)$mean
    len <- length(data.set(cpt))

    binList <- cl.iter(pts, len)

    if (length(binList) != length(ms)){
        ms <- lapply(binList, function(xx){ return( mean(data.set(cpt)[unlist(xx)]) ) })
    }

    lens <- sapply(binList, length)
    big <- ms[which(lens == max(lens))]
    FCs <- unlist(ms)/as.numeric(big)
    bad <- binList[which(FCs > k)]
    unlist(bad)
}

bin.unlist <- function(bins){
   cls <- lapply(bins, is.list)
   bins2 <- list()
   for (ii in 1:length(cls)){
       if (cls[[ii]]){
           vv <- lapply(bins[[ii]], unlist)
           bins2 <- c(bins2, vv)
       }
       else { bins2 <- c(bins2, bins[ii]) }
   }
   bins2
}

posList <- function(bins){
    ids <- sapply(bins, function(xx){
        if (length(xx) == 1 && xx == -1){ return(FALSE) }
        else { return(TRUE) }
    })
    return(bins[ids])
}

cl.iter <- function(vectr, max.x){
   bins <- lapply(vectr, cl, vectr=vectr, max.x=max.x)
#   browser()
   bins <- posList(bin.unlist(bins))
   bins
}
 
cl <- function(x, vectr, max.x){
  id <- which(x == vectr)
  if (x == 1){
    return(x)
  }
  ## if the first CPT is not 1, and the 2nd CPT is CPT.1 + 1, return 2 entries.
  if (x == vectr[1] & vectr[1] != 1 & length(vectr) != 1){
    if ((vectr[2] - x) == 1){ return(c(1:(x - 1),x)) }
    else { return(1:x) }
  }
  if (x == vectr[length(vectr)]){
    if (x != max.x){ return((x+1):max.x) }
    else if (x == max.x & length(vectr) > 1){ return((vectr[id-1] + 1):x) }
    else if (x == max.x & length(vectr) == 1){ return(c(unlist(1:(x-1)), unlist(x))) }
    else { return(x) }
  }
  else{
    if (((vectr[id+1] - x) == 1) & ((vectr[id+1] - vectr[id-1]) ==  2)){ return(x) }
    else if (((vectr[id+1] - x) == 1) & ((vectr[id+1] - vectr[id-1]) >  2)){
      return(-1) }
    else if (((vectr[id+1] - x) > 1) && ((vectr[id+1] - vectr[id-1]) >  2) &&
             ((id - 2) > 0) && ((x - vectr[id-1]) > 1)  &&
             ((x - vectr[id-2]) >  2)){
        return((x+1):(vectr[id+1])) }
    else { return(x:(vectr[id+1])) }
  }
}


lp <- function(dF, p=NULL){
  xnorms <- apply(dF, 2, function(x){
    pos <- x[which(x > 0)]
    if (is.null(p)){
      p = length(pos)
    }
    lnorm <- sum(abs(pos)^p)^(1/p)
  })
  return(xnorms)
}


