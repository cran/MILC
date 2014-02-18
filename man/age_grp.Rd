\name{age_grp}
\alias{age_grp}
\title{Define age group}
\description{
Returns the 5-year age group ("<40", "40-44", ..., or "85+" years) to which an individual belongs.
}
\usage{
age_grp(age)
}
\arguments{
  \item{age}{ continuous variable in years.
}
}
\value{ Factor indicating the 5-year age group.
}
\author{
Stavroula A. Chrysanthopoulou
}
\seealso{\code{\link{d_grp}}
}
\examples{
# 5-year age group of an individual 21 years old:
  age_grp(21)
# 5-year age group of an individual 77.5 years old:
  age_grp(77.5)
}
\keyword{Functions}