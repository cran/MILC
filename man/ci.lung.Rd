\name{ci.lung}
\alias{ci.lung}
\docType{data}
\title{ci.lung dataset
}
\description{
An R object of class "list" containing Cumulative Incidence Function estimates for lung cancer mortality, by age (5-years age group), gender ("male" or "female") and tumor stage at diagnosis (i.e., "localized", "regional", or "distant").
}
\usage{ci.lung}
\format{
The \emph{`ci.lung'} list comprises 68 components in total. Each component itself is an R object of class "list" with three components, i.e.:

[[1]] a vector with time points (in years)

[[2]] a vector with CIF estimates 

[[3]] a vector with variances of the CIF estimates
}
\source{
CIF estimates are derived using data obtained from the Surveillance, Epidemiology, and End Results (SEER) database covering the years 1973-2008.
}
\seealso{ \code{\link{current.other}, \link{former.other}, \link{never.other}}
}
\examples{
data(ci.lung)
ci.lung[1]
ci.lung[[1]][1]		# time (years)
ci.lung[[1]][2]		# CIF estimates
ci.lung[[1]][3]		# variances of CIF estimates
}
\keyword{Datasets}