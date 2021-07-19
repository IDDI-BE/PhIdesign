
#' @title Get BOIN decision rules for phase I design
#' @description Get BOIN decision rules for phase I design
#' @param phi target toxicity rate
#' @param phi1 the highest DLT rate that is deemed subtherapeutic (i.e., underdosing), 
#'     such that dose escalation should be made
#' @param phi2 the lowest DLT rate that is deemed overly toxic (i.e., overdosing), 
#'     such that dose de-escalation is required
#' @return a data.frame with elements
#' \itemize{
#' \item phi
#' \item phi1
#' \item phi2
#' \item lambda_e: if observed DLT rate<= lambda_e --> dose escalation
#' \item lambda_d: if observed DLT rate>= lambda_d --> dose de-escalation
#' }
#' @references Liu S, Yuan Y. Bayesian Optimal Interval Designs for Phase I Clinical Trials. 
#'     J R Stat Soc Ser C Appl Stat. 2015;64:507–23
#' @export
#' @examples
#' rule<-get_BOIN_rules(phi=0.3,phi1=0.6*0.3,phi2=1.4*0.3)

get_BOIN_rules<-function(phi,phi1,phi2){
  
  lambda_e<- log((1-phi1)/(1-phi ))/log((phi *(1-phi1))/(phi1*(1-phi )))
  lambda_d<- log((1-phi )/(1-phi2))/log((phi2*(1-phi ))/(phi *(1-phi2)))
  
  return(data.frame(phi=phi,phi1=phi1,phi2=phi2,lambda_e=lambda_e,lambda_d=lambda_d))

}