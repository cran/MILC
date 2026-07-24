\name{d_grp}
\alias{d_grp}
\title{Define smoking intensity group
}
\description{
Returns the smoking intensity group ("1-10", "11-20", or "21+" cigarettes) of an individual based on the average number of smoked per day.
}
\usage{
d_grp(d)
}
\arguments{
  \item{d}{ continuous variable indicating smoking intensity (number of cigarettes).
}
}
\value{ Factor indicating the smoking intensity group.
}
\author{
Stavroula A. Chrysanthopoulou}
\seealso{\code{\link{age_grp}}
}
\examples{
# Specify the smoking intensity group of an individual who smokes: 

# 1 cigarette/day:
  d_grp(1)

# 22 cigarettes/day:
  d_grp(22)
}
\keyword{Functions}