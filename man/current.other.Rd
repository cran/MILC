\name{current.other}
\alias{current.other}
\docType{data}
\title{current.other dataset
}
\description{
An R object of class "list" containing Cumulative Incidence Function estimates for other cause (not lung cancer) mortality, for current smokers by age (5-years age group), gender ("male" or "female") and smoking intensity group ("1-10", "11-20", or "21+" cigerettes) based on the average number of cigarettes smoked per day.
}
\usage{current.other}
\format{
The \emph{`current.other'} list comprises 66 components in total. Each component itself is an R object of class "list" with three components, i.e.:

[[1]] a vector with time points (in years)

[[2]] a vector with CIF estimates 

[[3]] a vector with variances of the CIF estimates

}
\source{
CIF estimates are derived using data obtained from the National Health Interview Survey (NHIS) database covering the years 1986-2004.
}
\seealso{ \code{\link{ci.lung}, \link{former.other}, \link{never.other}}
}
\examples{
data(current.other)
current.other[1]
current.other[[1]][1]	# time (years)
current.other[[1]][2]	# CIF estimates
current.other[[1]][3]	# variances of CIF estimates
}
\keyword{Datasets}
