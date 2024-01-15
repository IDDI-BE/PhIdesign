
CJ.dt <- function(X, Y) {
  .SD<-NULL
  k <- NULL
  X <- X[, c(k = 1, .SD)]
  data.table::setkey(X, k)
  Y <- Y[, c(k = 1, .SD)]
  data.table::setkey(Y, NULL)
  X[Y, allow.cartesian = TRUE][, `:=`(k, NULL)]
} 

#' @title Get Keyboard decision rules for phase I design
#' @description Get Keyboard decision rules for phase I design
#' @param N sample size
#' @param phi target toxicity rate
#' @param halfkey half length of key in keyboard design
#' @return a data.frame with elements
#' \itemize{
#' \item n: number of patients
#' \item x: observed DLTs
#' \item decision: "E" (escalate),"R" (retain) or "D" (de-escalate)
#' }
#' @details
#' A uniform Beta(1,1) prior is used to calculate the posterior Bayesian probability for each 'key'
#' @references Yan F, Mandrekar SJ, Yuan Y. Keyboard: A Novel Bayesian Toxicity Probability Interval Design for Phase I Clinical Trials. Clin Cancer Res. 2017; 23: 3994–4003
#' @export
#' @examples
#' Keyboard_rules <- get_Keyboard_rules(N=9,phi=0.3,halfkey=0.05)

get_Keyboard_rules<-function(N,phi,halfkey){
  
  n<-x<-key_low<-key_high<-keyprob<-decision<-NULL
    
  # Create dataset with all keys for each combination of "n" (number of patients) and "x" (number of successes)
  key <- data.table::as.data.table(data.frame(key=c(rev(seq(phi-halfkey,0,-2*halfkey)),seq(phi+halfkey,1,2*halfkey))))
  res <- lapply(1:N, FUN=function(n) cbind(n = n, x = 0:n)) # create dataset with for each N, all possible outcomes (nr of DLTs)
  res <- data.table::as.data.table(do.call(rbind,  res))
  
  res <- CJ.dt(res,key)
  res <- res[order(n,x,key)]
  
  res <- cbind(res[1:(dim(res)[1]-1),],res[2:(dim(res)[1]),"key"])
  names(res)<-c("n","x","key_low","key_high")
  res <- res[key_low<key_high] 
  
  res <- res[,keyprob:=pbeta(key_high,shape1=x+1,shape2=n-x+1)-pbeta(key_low,shape1=x+1,shape2=n-x+1)]
  res <- res[,sum:=sum(keyprob),by=list(n,x)]
  res <- res[,max:=max(keyprob),by=list(n,x)]
  
  Decision <- res[keyprob==max]
  
  Decision<-Decision[,decision:=data.table::fifelse(key_high <(phi+halfkey),"E",
                                                    data.table::fifelse(key_high >(phi+halfkey),"D",
                                                                        data.table::fifelse(key_high==(phi+halfkey),"R",NA_character_)))]
  Decision<-Decision[,c("n","x","decision")]
  
  return(Decision)
}