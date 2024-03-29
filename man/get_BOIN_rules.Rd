% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_BOIN_rules.R
\name{get_BOIN_rules}
\alias{get_BOIN_rules}
\title{Get BOIN decision rules for phase I design}
\usage{
get_BOIN_rules(phi, phi1, phi2)
}
\arguments{
\item{phi}{target toxicity rate}

\item{phi1}{the highest DLT rate that is deemed subtherapeutic (i.e., underdosing), 
such that dose escalation should be made}

\item{phi2}{the lowest DLT rate that is deemed overly toxic (i.e., overdosing), 
such that dose de-escalation is required}
}
\value{
a data.frame with elements
\itemize{
\item phi
\item phi1
\item phi2
\item lambda_e: if observed DLT rate<= lambda_e --> dose escalation
\item lambda_d: if observed DLT rate>= lambda_d --> dose de-escalation
}
}
\description{
Get BOIN decision rules for phase I design
}
\examples{
rule<-get_BOIN_rules(phi=0.3,phi1=0.6*0.3,phi2=1.4*0.3)
}
\references{
Liu S, Yuan Y. Bayesian Optimal Interval Designs for Phase I Clinical Trials. 
    J R Stat Soc Ser C Appl Stat. 2015;64:507–23
}
