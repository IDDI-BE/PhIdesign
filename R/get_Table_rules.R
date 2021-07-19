
#' @title Get table with summary of decision rules for phase I design
#' @description Get table with summary of decision rules for phase I design (BOIN or Keyboard design)
#' @param N sample size
#' @param phi target toxicity rate
#' @param phi1 the highest DLT rate that is deemed subtherapeutic (i.e., underdosing), 
#'     such that dose escalation should be made
#' @param phi2 the lowest DLT rate that is deemed overly toxic (i.e., overdosing), 
#'     such that dose de-escalation is required
#' @param maxtox P(DLT rate>phi|c,n)<=maxtox as value for maxprob in \code{\link{get_Elim_rules}}
#' @param halfkey half length of key in keyboard design
#' @param design "BOIN" or "Keyboard"
#' @return a data.frame with elements
#' \itemize{
#' \item n: number of patients
#' \item E: max(observed DLTs) at n --> decision to escalate dose
#' \item R: range(observed DLTs) at n --> decision to retain dose
#' \item R: range(observed DLTs) at n --> decision to retain dose
#' \item lambda_e: if observed DLT rate<= lambda_e --> dose escalation
#' \item lambda_d: if observed DLT rate>= lambda_d --> dose de-escalation
#' }
#' @references Liu S, Yuan Y. Bayesian Optimal Interval Designs for Phase I Clinical Trials. 
#'     J R Stat Soc Ser C Appl Stat. 2015;64:507–23
#' @export
#' @examples
#' BOIN_table<-get_Table_rules(N=21,phi=0.3,phi1=0.6*0.3,phi2=1.4*0.3,maxtox=0.95,design="BOIN")
#' KEYB_table<-get_Table_rules(N=21,phi=0.3,maxtox=0.95,halfkey=0.05,design="Keyboard")

get_Table_rules<-function(N,phi,phi1=NULL,phi2=NULL,maxtox,halfkey=NULL,design){

  E_max<-x<-n<-E<-D_min<-D<-R_min<-R<-R_max<-decision<-NULL
  
  if (design=="BOIN"){
  
    BOIN<-get_BOIN_rules(phi=phi,phi1=phi1,phi2=phi2)
    lambda_e<-BOIN$lambda_e
    lambda_d<-BOIN$lambda_d
    
    res   <- lapply(1:N, FUN=function(n) cbind(n = n, x = 0:n)) # create dataset with for each N, all possible outcomes (nr of DLTs)
    res   <- data.table::as.data.table(do.call(rbind, res))
    res$p <- res$x/res$n
    res$E <- ifelse(res$p<=lambda_e,1,0)                 # Escalate
    res$D <- ifelse(res$p>=lambda_d,1,0)                 # De-escalate
    res$R <- ifelse(lambda_e<res$p & res$p<lambda_d,1,0) # Retain
    
    res <- res[, E_max:=max(x), by = list(n,E)] # Get max(x) to escalate
    res <- res[, D_min:=min(x), by = list(n,D)] # Get min(x) to de-escalate
    res <- res[, R_min:=min(x), by = list(n,R)] # Get min(x) to retain
    res <- res[, R_max:=max(x), by = list(n,R)] # Get max(x) to retain
    
    Escalate      <- unique(res[E==1,c("n","E_max")],by="n")
    DeEscalate    <- unique(res[D==1,c("n","D_min")],by="n")
    Retain        <- unique(res[R==1,c("n","R_min","R_max")],by="n")

  }
  
  else if (design=="Keyboard"){
    
    res <- get_Keyboard_rules(N=N,phi=phi,halfkey=halfkey)
    
    res <- res[, E_max:=max(x), by = list(n,decision)] # Get max(x) to escalate
    res <- res[, D_min:=min(x), by = list(n,decision)] # Get min(x) to de-escalate
    res <- res[, R_min:=min(x), by = list(n,decision)] # Get min(x) to retain
    res <- res[, R_max:=max(x), by = list(n,decision)] # Get max(x) to retain
    
    Escalate      <- unique(res[decision=="E",c("n","E_max")],by="n")
    DeEscalate    <- unique(res[decision=="D",c("n","D_min")],by="n")
    Retain        <- unique(res[decision=="R",c("n","R_min","R_max")],by="n")
  }
  
  Retain$Retain <-ifelse(Retain$R_min==Retain$R_max,Retain$R_min,
                         ifelse(Retain$R_min!=Retain$R_max,paste(Retain$R_min,Retain$R_max,sep="-"),NA))
  alln          <- data.table::data.table(n=1:N)
  Decision_ret  <- merge(alln,Retain[,c("n","Retain")],all=T); names(Decision_ret)[names(Decision_ret)=="Retain"]<-"R"
  Decision      <- merge(Escalate,Decision_ret,by="n"); names(Decision)[names(Decision)=="E_max"]<-"E"
  Decision$D    <- DeEscalate$D_min
  Decision$Elim <- get_Elim_rules(Nmax=N,phi=phi,maxprob=maxtox)$c_STOP
  
  return(Decision)
}