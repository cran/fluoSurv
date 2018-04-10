#' fluoSurv: A package for estimating insect survival data from spectrophotometry measurements
#'
#' The fluoSurv package provides functions to import fluorescence data as exported from a BioTek
#' microplate reader and functions to estimate insect survival from these fluorescence data.
#'
#' @section The function to import data:
#' \link{read.kinetic}
#'
#' @section The function to estimate survival:
#' \link{estimate.LT}
#'
#' @docType package
#' @name fluoSurv
NULL

#' Extract data for a single well
#'
#' @param data A \code{data.frame} containing fluorescence measurements for the 96 wells
#' of a plate and several time values. See \cite{read.kinetics} for a way to obtain such
#' a \code{data.frame} from data files produced by a Biotech plate-reader.
#'
#' @details The number of measurements might differ between wavelengths. In the Biotek reader used here, this can happen if fluorescence value exceeds the maximum value of 10^6. \code{NA} are then added to the \code{data.frame}.
#'
#' @return A \code{data.frame}
#' @export
#'
extract.well <- function (data) {
  if(! ("ID_read" %in% colnames(data)) ) return(NULL)
  ## ID of the different reads (i.e. wave-length / gain combination)
  liste.ID.read <- levels(data$ID_read)
  dtmp3 <- NULL
  for(i in 1:length(liste.ID.read)) {
    cat(liste.ID.read[i],"\n")
    dtmp2 <- subset(data,data$ID_read==liste.ID.read[i])

    if(is.null(dtmp3)) dtmp2 <- dtmp2[,c("num","t","value")]

    else dtmp2 <- dtmp2[,c("t","value")]
    colnames(dtmp2)[ncol(dtmp2)-1:0] <- paste(c("t","value"),liste.ID.read[i],sep="_")
    if(is.null(dtmp3)) dtmp3 <- dtmp2
    else {
      if(nrow(dtmp2)==nrow(dtmp3)) {
        dtmp3 <- cbind(dtmp3,dtmp2)
      } else {
        if(nrow(dtmp2)<nrow(dtmp3)) {
          dtmp4 <- matrix(rep(c(NA,NA),nrow(dtmp3)-nrow(dtmp2)),ncol=2)
          dtmp4 <- as.data.frame(dtmp4)
          colnames(dtmp4) <- colnames(dtmp2)
          dtmp3 <- cbind(dtmp3,rbind(dtmp2,dtmp4))
          cat("WARNING! NA are added to ",liste.ID.read[i],"because it contains",nrow(dtmp2),"measurement where ",nrow(dtmp3),"are expected!\n")
        } else {
          cat("WARNING!",liste.ID.read[i],"is ignored because it contains",nrow(dtmp2),"measurement where ",nrow(dtmp3),"are expected!\n")
        }
      }
    }
  }
  return(dtmp3)
}

###########################################
##### GLM threshold model  ################
###########################################

## moving average function
ma <- function(x,n=50){filter(x,rep(1/n,n), sides=2)}


#' Computes when fluorescence exceeds a given threshold value
#' @param t The time value
#' @param x The fluorescence value
#' @param min.t The time value after which threshold value is searched
#' @param threshold Threshold value, as a proportion above the maximum intensity value observed before min.t
#' @param n Width of the moving average window used to smooth signal
#' @examples
#'
#' data(galleria)
#' d <- subset(galleria,!is.na(value))
#' l  <- lapply(split(d,d$well),extract.well)   #complete kinetics for each well
#'
#' with(l[["A3"]],plot(t_2_485_535,log(value_2_485_535,10),type="l"))
#' with(l[["A3"]],abline(v=when.threshold(t_2_485_535,value_2_485_535)))
#'
#' @export
when.threshold <- function(t,x,min.t=5,threshold=0.1,n=50){
  if(is.null(x) || is.null(t)) return(NA)
  if(!any(t>min.t)) return(NA)

  ma.x <- ma(x,n=n)
  l <- ma.x > max(ma.x[t<=min.t],na.rm=T)*(1+threshold)
  if(!any(l[!is.na(l)])) return(NA)
  z <- min(t[which(l)],na.rm=T)
  if(z<=min.t) return(NA)
  return(z)
}


#' Estimate variance in signal
#'
#' @param l A list object with each element being measurement for a sample
#' @param var  A character string corresponding to the name of the signal used to estimate lethal time
#' @param min.t The time after which signal will be used to compute offset
#' @param max.t The time before which signal will be used to compute offset
#'
#' @return A value that can be used as an offset in \cite{estimate.LT}
#' @export
estimate.offset <- function(l,var,min.t=NULL,max.t=NULL) {
  var.t <- paste("t",var,sep="_")
  var.y <- paste("value",var,sep="_")
  if(!("list" %in% class(l))) l <- list(l)
  v.sd <- sapply(l,function(dtmp) {
    if(!(var.t %in% colnames(dtmp))) {
      print(head(dtmp))
      stop(paste(var.t,"not found in data!"))
    }
    if(!(var.y %in% colnames(dtmp))) {
      print(head(dtmp))
      stop(paste(var.y,"not found in data!"))
    }

    t  <- dtmp[,var.t] ## time values
    y  <- dtmp[,var.y] ## signal values
    fl <- !is.na(y)
    if(!is.null(min.t)) fl <- fl & t>= min.t
    if(!is.null(max.t)) fl <- fl & t<= max.t
    mean(diff(log(y[fl],10))^2)
  })

  offset <- 1/median(v.sd)
  return(offset)
}

# Returns the summed likelihoods of two Gamma glm where y^2 (with y the differentiated signal)
#is supposed to be constant. One model is for t<=td, the second for t>td. Estimates of y^2 can be
#provided for both models as offsets. If they are not, variance is estimated (as the intercept of the model).
fit.threshold.model <- function(td,y,t,offset.alive=NULL,offset.dead=NULL) {
  loglik <- NA

  l <- !is.na(y)  & t<=td
  n <- sum(l, na.rm=TRUE)
  if(n>1) {
    # Model 1 for t<=td
    if(is.null(offset.alive)) m1 <- try(stats::glm(y[l]~1,family=stats::Gamma())) # variance estimate is not provided
    else m1 <- try(glm(y[l]~ -1,offset=rep(offset.alive,n),family=Gamma()))       # variance estimate is provided and used as offset
    if("try-error" %in% class(m1)) return(NA)
    loglik <- logLik(m1)
  }

  l <- !is.na(y)  & t>td
  n <- sum(l, na.rm=TRUE)
  if(n>1) {
    # Model 2 for t>td
    if(is.null(offset.dead)) m2 <- try(glm(y[l]~1,family=Gamma()))                       # variance estimate is not provided
    else m2 <- try(stats::glm(y[l]~ -1,offset=rep(offset.dead,n),family=stats::Gamma())) # variance estimate is provided and used as offset
    if("try-error" %in% class(m2)) return(NA)
    if(!is.na(loglik)) loglik <- loglik + stats::logLik(m2)
    else loglik <- stats::logLik(m2)
  }

  return(loglik)
}

#' Estimation of time to death
#'
#' @param y The signal to analize
#' @param t The time values
#' @param threshold.value Detection threshold below which signal is considered as pure noise
#' @param offset.dead  Offset derived from variance in signal in dead insects
#' @param offset.alive Offset derived from variance in signal in living insects
#' @param verbose If true additionnal informations on computations are displayed
#' @param ndeps Step size used in optim to estimate parameters
#'
#' @details The model adjusted in this procedure assumes that random variance in signal drops after insect death,
#'          because insect has ceased to move. The time value at which this drop occurs can therefore be used as
#'          an estimate of lethal time. First guess of variance estimates can be provided as offset.dead (for dead
#'          insects) or offset.alive (for living insects).
#'
#'
#' @return A vector with TL estimate for the given sample, corresponding log-likelihood and number of values used in the computation.
#' @importFrom utils write.table head
#' @importFrom stats glm Gamma coef filter logLik median
#' @export
#'
#' @examples
#'
#' ##Loading data
#' data("galleria")
#' ## dataset may contain NA if microplate reader has been stoped before the programmed
#' ## end of the experiment
#' d <- subset(galleria,!is.na(value))
#' l  <- lapply(split(d,d$well),extract.well)
#'
#' data(setup)
#' setup <- setup[match(setup$well,names(l)),]
#'
#' ## Computes rough estimates of variance in autofluorescence signal for dead and living insects
#' ## These values serve as initial guess to fit the model.
#' offset.alive <- estimate.offset(l,"1_330_405",min.t=1,max.t=5)
#'     # all insects are assumed to be alive during the
#'     # first five hours that follow injection
#'
#' offset.dead <- estimate.offset(l[which(setup$dead==1)],"1_330_405",min.t=72-5)
#'     # insects that were dead at the end of the experiment are
#'     # assumed to have died earlier than 5 hours before then
#'     # end of the experiment
#'
#'
#' ## LT estimation or for a single well
#' ## Check out well D9, to see what happens when an insect stayed alive.
#' well <- "A3"
#' plot(log(value_1_485_535,10)~t_1_485_535,type="l",col="green",ylim=c(2,5),data=l[[well]])
#' points(log(value_1_330_405,10)~t_1_330_405,type="l",col="gray",data=l[[well]])
#'     # Rough estimate obtained using no offsets
#' est1 <- with(l[[well]], estimate.LT(value_1_330_405,t_1_330_405,threshold.value=3))
#' abline(v=est1[["LT"]],lty=2,col="red")
#'     # Much better estimate obtained using offset for dead insects
#' est2 <- with(l[[well]], estimate.LT(value_1_330_405,t_1_330_405,
#'                                          offset.dead = offset.dead,threshold.value=3))
#' abline(v=est2["LT"],lty=3,col="red")
#'     # Using offset.alive does not change anything to the estimate for well A3
#'     # It may help for insect that have a larger variance in signal than others even after death
#' est3 <- with(l[[well]], estimate.LT(value_1_330_405,t_1_330_405,
#'                                          offset.dead = offset.dead,
#'                                          offset.alive = offset.alive,threshold.value=3))
#' abline(v=est3["LT"],col="red")
#'     # Detection of significant GFP fluorescence (i.e. log fluorescence exceed
#'     # by 5% the maximum value observed during the first five hours)
#' with(l[[well]],abline(v=when.threshold(t_1_485_535,log(value_1_485_535,10),
#'                                                             threshold=0.1),col="green"))
#'
#'
#' ##LT estimation for all wells
#' if(FALSE) { #example takes time! Set to TRUE if you want to run it
#'    res <- sapply(l,function(x) estimate.LT(x$value_1_330_405,x$t_1_330_405,
#'                                                threshold.value=2,offset.dead=offset.dead,
#'                                                offset.alive=offset.alive))
#'   res <- as.data.frame(t(res))
#' ##Adds LT estimates to the experimental setup data.frame
#'   setup <- cbind(setup,res[match(setup$well,rownames(res)),])
#' ## Time of injection is added to LT, so that LT really corresponds to time post injection
#'   time <- with(setup,strptime(as.character(time_injection),format="%H:%M:%S"))
#'   time <- as.numeric((max(time)-time)/(60^2))
#'   setup$LT <- setup$LT+time
#'
#'
#' ## Survival curves by dilution of bacterial culture injected
#'   library(survival)
#'   plot(survfit(Surv(LT,dead)~dilution,data=setup),
#'      lwd=c(3:1,1),lty=c(1,1,1,2),
#'      xlab="hours post injection",ylab="proportion of surviving insects")
#'   abline(h=0.5,col="red")
#'   legend("topright",lwd=c(1,3:1),lty=c(2,1,1,1),legend=c("LB",10^(3:1)))
#'
#' ## When does scepticemia start?
#'   res <- sapply(l,function(x) when.threshold(x$t_1_485_535,log(x$value_1_485_535,10),
#'                                                                     threshold=0.1))
#'   setup$T_gfp <- res[match(setup$well,names(res))]
#'   setup$T_gfp <- setup$T_gfp + time
#'
#' ## Relation between time of death and moment when scepticemia is detected.
#' ## Only points where scepticemia has been detected are represented here.
#'   plot(LT~T_gfp,data=setup,col=ifelse(dead,1,2),pch=as.numeric(dilution))
#'   abline(0,1)
#'   with(setup,legend("topleft",legend=levels(dilution),pch=1:4))
#'   # Most insects have died after scepticemia has started.
#'   }
#'
estimate.LT <- function(y,t,threshold.value=NULL,offset.dead=NULL, offset.alive=NULL, verbose=F,ndeps=0.01) {
  yy <- diff(log(y,10))^2
  if(is.null(threshold.value)) l <- !is.na(yy) & yy>0
  else l <- !is.na(yy) & yy>0 & log(y[-1],10)>threshold.value
  yy <- yy[l]
  tt <- (t[-1])[l]
  n <- length(tt)

  res <- c(rep(NA,2),n)
  names(res) <- c("LT","logLik","n")

  if(is.null(offset.dead)) {
    sol <- try(stats::optim(function(x) -fit.threshold.model(x,yy,tt),par=tt[n-4],method="L-BFGS-B",lower=tt[1],upper=tt[n],control=list(ndeps=ndeps)))
    if(verbose)  print(sol)
    if(!("try-error" %in% class(sol))){
      res[1] <- sol$par
      res[2] <- sol$value
    }
    return(res)
  }

  #first estimation of td (fixed variance after death)
  sol1 <- try(stats::optim(function(x) -fit.threshold.model(x,yy,tt,offset.dead=offset.dead,offset.alive=offset.alive),par=tt[n-40],method="L-BFGS-B",
                    lower=tt[1],upper=tt[n],control=list(ndeps=ndeps)),silent=T)
  if(verbose) print(sol1)
  if("try-error" %in% class(sol1)) return(res)

  #re-estimation of variance after death
  m <- glm(yy~ I(tt>sol1$par),family=Gamma())
  if(verbose) cat(logLik(m),"\n")

  #re-estimation of td with new variance estimate
  sol2 <- try(stats::optimize(function(x) -fit.threshold.model(x,yy,tt,offset.alive=coef(m)[1],offset.dead=sum(coef(m))),lower=tt[3],upper=tt[n-2]))
  if(verbose) print(sol2)
  if(!("try-error" %in% class(sol2))) {
    res[1:2] <- c(sol2$minimum,-sol2$objective)
  }

  return(res)
}
