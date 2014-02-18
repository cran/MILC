\name{former.other}
\alias{former.other}
\docType{data}
\title{former.other dataset
}
\description{
An R object of class "list" containing Cumulative Incidence Function estimates for other cause (not lung cancer) mortality, for former smokers by age (5-years age group), gender ("male" or "female") and smoking intensity group ("1-10", "11-20", or "21+" cigarettes) based on the average number of cigarettes smoked per day.
}
\usage{former.other}
\format{
The \emph{`former.other'} list comprises 22 components in total. Each component itself is an R object of class "list" with three components, i.e.:

[[1]] a vector with time points (in years)

[[2]] a vector with CIF estimates 

[[3]] a vector with variances of the CIF estimates

}
\source{
CIF estimates are derived using data obtained from the National Health Interview Survey (NHIS) database covering the years 1986-2004.
}
\seealso{ \code{\link{ci.lung}, \link{current.other}, \link{never.other}}
}
\examples{
data(former.other)
former.other[1]
former.other[[1]][1]	# time (years)
former.other[[1]][2]	# CIF estimates
former.other[[1]][3]	# variances of CIF estimates
}
\keyword{Datasets}