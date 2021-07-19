
#' @title Get dose elimination rules for phase I design
#' @description Get Dose elimination rules for phase I design, based on Bayesian posterior probability
#' @param Nmax maximum number of patients
#' @param phi target DLT rate
#' @param maxprob P(DLT rate>phi|c,n)<=maxprob
#' @return a data.frame with elements
#' \itemize{
#' \item phi: target toxicity rate
#' \item n: number of patients at which toxicity is evaluated
#' \item c: max number of DLTs at which P(DLT rate>phi)<=maxprob
#' \item c_STOP: c+1, i.e. min number of DLTs at which P(DLT rate>phi|c_STOP,n)>maxprob
#' \item prob_Bayes: P(DLT rate>phi|c,n)
#' \item prob_Obs: c/n or the allowed observed toxicity rate with no stopping
#' }
#' @details
#' A uniform Beta(1,1) prior is used to calculate the posterior Bayesian probability
#' @export
#' @examples
#' rule<-get_Elim_rules(Nmax=21,phi=0.3,maxprob=0.95)


get_Elim_rules<-function(Nmax,phi,maxprob){
  
  result<-data.frame(phi=phi,n=1:Nmax,c=rep(NA,Nmax),c_STOP=rep(NA,Nmax),prob_Bayes=rep(NA,Nmax),prob_Obs=rep(NA,Nmax))
  
  for (i in 1:Nmax){
    c        <- i+1
    postprob <- 1
    while (postprob>maxprob){
      c<-c-1
      postprob<-1-pbeta(phi,shape1=c+1,shape2=i-c+1)
    }
    result[i,"c"]         <-c
    result[i,"prob_Bayes"]<-round(postprob,2)
  }
  
  result$prob_Obs <- round(result$c/result$n,2)
  result$c_STOP   <- (result$c)+1
  return(result)
  
}