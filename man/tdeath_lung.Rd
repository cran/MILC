\name{tdeath_lung}
\alias{tdeath_lung}
\title{Predict the age at death from lung cancer
}
\description{Function to predict the age (years) at which a person may die of lung cancer given gender, age and tumor stage at diagnosis.
}
\usage{
tdeath_lung(u1, u2, covs_lung)
}
\arguments{
  \item{u1, u2}{random numbers from Unif[0,1] required for the simulation
}
  \item{covs_lung}{3-dimensional vector with values for the covariates related to death from lung cancer, i.e., tumor stage at diagnosis ("localized", "regional", or "distant"), age (years) at diagnosis, and gender ("male" or "female")} 
}
\value{An R-object of class "list" with the following six components:

[[1]]: random number u1 used in the simulation

[[2]]: random number u2 used in the simulation 

[[3]]: index number of the time interval

[[4]]: age (years) at death from lung cancer

[[5]]: whether the person died of lung cancer("Yes") or not("No")

[[6]]: R-object of class "list" with the relevant CIF estimates
}
\author{
Stavroula A. Chrysanthopoulou}
\note{
Components [[1]]-[[3]] and [[6]] are returned for testing purposes only.
}
\seealso{\code{\link{ci.lung}, \link{tdeath_other}}
}
\examples{
# Predict the age at death from lung cancer for a woman, diagnosed with 
# "localized" tumor at the age of 89 years
data(ci.lung)
d.lung <- tdeath_lung(runif(1), runif(1), c("localized", 89, "female"))
d.lung[[1]]
d.lung[[2]]
d.lung[[3]]
d.lung[[4]]
d.lung[[5]]
d.lung[[6]]
}
\keyword{Functions}