#' Adjusted p value for testing H_J using Dunnett method
#'
#'This functions provides the adjusted p value for testing a family-wise hypothesis H_J based on the p values of individual raw p values in J.
#'
#' @param p A vector of individual raw p values in J.
#' @param cr correlation structure.
#'
#' @return Adjusted p values using Dunnett procedure
#'
#' @importFrom mvtnorm pmvnorm
#' @importFrom stats qnorm
#' @export
#'
#' @examples
#' p = c(0.01, 0.02, 0.03, 0.013)
#' cr <- matrix(0.5, 4, 4); diag(cr) <- 1
#' dunnett(p, cr)
dunnett <- function(p, cr){
  num_trt <- length(p)
  # equal weight
  w <- rep(1, num_trt)/num_trt
  padjusted <- p.dunnet.ph23(p = p, cr = cr,w = w, upscale = FALSE)
  ans <- min(min(padjusted), 1)
  return(ans)
}

## direct copy of gMCP:::p.dunnet, but added seed = 20240716 for pmvnorm so that the results are reproducible.
p.dunnet.ph23 <- function(p,cr,w,upscale, alternatives="less"){
  if(length(cr)>1){
    conn <- conn.comp.ph23(cr)
  } else {
    conn <- 1
  }
  twosided <- alternatives==rep("two.sided", length(w))
  lconn <- sapply(conn,length)
  conn <- lapply(conn,as.numeric)
  e <- sapply(1:length(p),function(i){
    sum(sapply(conn,function(edx){
      if(length(edx)>1){
        if (upscale=="o3") {
          return((1-mvtnorm::pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w))))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]*sum(w)))))),
            corr=cr[edx,edx], seed = 20240716, abseps=10^-5)))
        } else {
          return((1-mvtnorm::pmvnorm(
            lower=ifelse(twosided[edx],qnorm(pmin(1,(w[edx]*p[i]/(w[i])))/2),-Inf),
            upper=ifelse(twosided[edx],qnorm(1-pmin(1,(w[edx]*p[i]/(w[i])))/2),qnorm(1-pmin(1,(w[edx]*p[i]/(w[i]))))),
            corr=cr[edx,edx], seed = 20240716,abseps=10^-5))/ifelse(upscale,1,sum(w)))
        }
      } else {
        if(upscale=="o3" || !upscale){
          return((w[edx]*p[i]/(w[i]*sum(w))))
        } else {
          return((w[edx]*p[i]/(w[i])))
        }
      }
    }))})

  e <- pmin(e,1)
  e
}

## direct copy of gMCP:::conn.comp
conn.comp.ph23 <- function(m){
  N <- 1:ncol(m)
  M <- numeric(0)
  out <- list()
  while(length(N)>0){
    Q <- setdiff(N,M)[1]
    while(length(Q)>0){
      w <- Q[1]
      M <- c(M,w)
      Q <- setdiff(unique(c(Q,which(!is.na(m[w,])))),M)
    }
    out <- c(out,list(M))
    N <- setdiff(N,M)
    M <- numeric(0)
  }
  return(out)
}
