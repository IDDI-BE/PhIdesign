
#' @title Get 'modified Finonacci' sequence for Phase I design
#' @description Get 'modified Finonacci' sequence for Phase I design
#'     Other decreasing step schemes can be specified
#' @param Ndose number of dose levels
#' @param first_dose First dose level
#' @param fibonacci decreasing step scheme, default is c(1,2,1.67,1.5,rep(1.33,Ndose-4)
#' @return a data.frame with elements
#' \itemize{
#' \item dose_level: dose level
#' \item dose: actual dose to be used
#' \item incr_ratio: The ratio of successive doses
#' }
#' @references Penel N, Kramar A. What does a modified-Fibonacci dose-escalation
#    actually correspond to? BMC Medical Research Methodology 2012, 12:103
#' @export
#' @examples
#' Fibonacci<-Get_Fibonacci(Ndose=6,first_dose=0.0017)

Get_Fibonacci<-function(Ndose,first_dose, fibonacci=c(1,2,1.67,1.5,rep(1.33,Ndose-4))){

  dose<-c(first_dose,rep(NA,Ndose-1))

  for (i in 2:Ndose){
    dose[i]<-dose[i-1]*fibonacci[i]
  }

  return(data.frame(dose_level=c(1:Ndose),dose=dose,incr_ratio=fibonacci))
}
