
#' @title Get i3+3 decision rules for phase I design
#' @description Get i3+3 decision rules for phase I design
#' @param N sample size
#' @param phi target toxicity rate
#' @param phi1 the highest DLT rate that is deemed subtherapeutic (i.e., underdosing), 
#'     such that dose escalation should be made
#' @param phi2 the lowest DLT rate that is deemed overly toxic (i.e., overdosing), 
#'     such that dose de-escalation is required
#' @return a data.frame with elements
#' \itemize{
#' \item n: number of patients
#' \item x: observed DLTs
#' \item decision: "E" (escalate),"R" (retain) or "D" (de-escalate)
#' }
#' @references Liu M et al. The i3+3 design for phase I clinical trials. 
#'     Journal of Biopharmaceutical Statistics 2020;30:294-304
#' @export
#' @examples
#' rule<-get_i33_rules(N=6,phi=0.3,phi1=0.25,phi2=0.35)

get_i33_rules<-function(N,phi,phi1,phi2){
  
  p1<-p2<-NULL
  
  # Create dataset with all combinations of "n" (number of patients) and "x" (number of successes)
  res <- lapply(1:N, FUN=function(n) cbind(n = n, x = 0:n)) # create dataset with for each N, all possible outcomes (nr of DLTs)
  res <- data.frame(do.call(rbind,  res))
  res$count<-1:dim(res)[1]
  
  res$p1=res$x/res$n
  res$p2<-ifelse(res$x!=0,(res$x-1)/res$n,NA)
  res$decision<-rep(NA,dim(res)[1])
  
  for (i in 1:dim(res)[1]){
    p1<-res[i,"p1"]
    p2<-res[i,"p2"]
    
    if (p1<phi1)             {res[i,"decision"]<-"E"} else
    if (p1>=phi1 & p1<=phi2) {res[i,"decision"]<-"R"} else
    if (!is.na(phi2)){
      if (p1>phi2){
        if (p2<phi1)             {res[i,"decision"]<-"R"} else
        if (p2>=phi1 & p2<=phi2) {res[i,"decision"]<-"D"} else
        if (p2>phi2)             {res[i,"decision"]<-"D"}
      }
    }
  }
  Decision<-res[,c("n","x","decision")]
  return(Decision)
}
