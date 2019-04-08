#' @importFrom stats smooth.spline predict
NULL

##' Glog fit
##'
##' @param x Values to convert.
##' @param a glog offset.
##' @param inverse To apply glog or reverse back to original values.
##' @export

glogFit <- function (x, a=1, inverse=FALSE) {
  if (inverse) {
    out <- 0.25 * exp(-x) * (4 * exp(2 * x) - (a * a))
  }else{
    out <- log((x + sqrt(x^2 + a^2))/2)
  }
  return(out)
}


##' Cubic smoothing spline fitting function for batch effect correction
##'
##' @param x QC sample measurement order indices.
##' @param y Signal intensities.
##' @param newX Range of sample indices to predict.
##' @param log TRUE of FALSE to perform the signal correction fit on the log scaled data. Default is TRUE.
##' @param a log scaling offset.
##' @param spar Smoothing parameter of the fitted curve. Should be in the range 0 to 1. If set to 0 it will be estimated using leave-one-out cross-validation.
##' @export

splineSmoother <- function(x, y, newX, log=log, a=1, spar) {
  if(log == TRUE) {
    y <- glogFit(y, a=a)
  }
  if (spar == 0) {
    sp.obj <- smooth.spline(x, y, cv=TRUE)
  } else {
    sp.obj <- smooth.spline(x, y, spar=spar)
  }

  valuePredict <- predict(sp.obj, newX)

  if(log==TRUE) {
    valuePredict$y <- glogFit(valuePredict$y, a=a, inverse=TRUE)
  }

  return(valuePredict$y)
}


##' Signal batch correction wrapper function
##'
##' @param id Feature id.
##' @param qcData Data frame with QC samples.
##' @param qcBatch QC batch numbers.
##' @param qcOrder QC measurement order.
##' @param spar Spline smoothing parameter. Should be in the range 0 to 1. If set to 0 it will be estimated using leave-one-out cross-validation.
##' @param log TRUE of FALSE to perform the signal correction fit on the log scaled data. Default is TRUE.
##' @param batch A vector indicating the batch each sample was measured in. If only one batch was measured then all values should be set to 1.
##' @param minQC Minimum number of QC samples required for signal correction.
##' @param order A numeric vector indicating the order in which samples were measured.
##' @export

sbcWrapper <- function(id, qcData, order, qcBatch, qcOrder, log=log, spar=spar, batch=batch, minQC) {
  out <- tryCatch ({

    # Measurment order can be non consecutive numbers as well
    maxOrder <- length(order)

    cat (id,"\n")
    newOrder <- 1:maxOrder
    subData <- qcData[id, ]

    nbatch <- unique(qcBatch)

    outl <- matrix(nrow=maxOrder, ncol=length(nbatch))

    for (nb in 1:length(nbatch)) {
      x <- qcOrder[qcBatch == nbatch[nb]]
      y <- subData[qcBatch == nbatch[nb]]
      NAhits <- which(is.na(y))

      if (length(NAhits) > 0) {
        x <- x[-c(NAhits)]
        y <- y[-c(NAhits)]
      }

      if (length(y) >= minQC) {
        outl[,nb] <- splineSmoother(x=x, y=y, newX=newOrder, log=log, a=1, spar=spar)
      } else {
        outl[,nb] <- rep(NA, maxOrder)
      }
    }

    outp <- rep(NA, nrow(outl))
    for (nb in 1:length(nbatch)) {
      range <- c(batch == nbatch[nb])
      outp[range] <- outl[range, nb]
    }

    outp
  },

  error = function(e){
    print (e, "\n")
    return(NULL)
  },

  warning = function(e){
    print (e, "\n")
    return(NULL)
  })

  return(outp)
}
