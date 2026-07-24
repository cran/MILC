\name{tdeath_other}
\alias{tdeath_other}
\title{Predict the age at death from a cause other than lung cancer
}
\description{Function to predict the age (years) at which a person may die from a cause other than lung cancer given age, gender and smoking intensity, when relevant.
}
\usage{
tdeath_other(u1, u2, status, covs_other)
}
\arguments{
  \item{u1, u2}{random numbers from Unif[0,1] required for the simulation
}
  \item{status}{smoking status ("never", "former", or "current" smoker)
}
  \item{covs_other}{3-dimensional vector with values for the covariates (other than smoking status) related to death from other causes, i.e., age (years) at the beginning of the prediction period, gender, smoking intensity expressed as average number of cigarettes smoked per day.
}
}
\value{
An R-object of class "list" with the following six components:

[[1]]: random number u1 used in the simulation

[[2]]: random number u2 used in the simulation 

[[3]]: index number of the time interval

[[4]]: time interval at which death from other causes may occur

[[5]]: age (years) at death from cause other than lung cancer

[[6]]: R-object of class "list" with the relevant CIF estimates
}
\note{
Components [[1]]-[[4]] and [[6]] are returned for testing purposes only.
}
\author{
Stavroula A. Chrysanthopoulou}
\seealso{\code{\link{current.other}, \link{former.other}, \link{never.other}, \link{tdeath_lung}}
}
\examples{
# Predict the age at death from a cause other than lung cancer for a man 52 years old, 
# who have never smoked.
data(current.other, former.other, never.other)
d.other <- tdeath_other(runif(1), runif(1), "never", c(52, "male", NA))
d.other[[1]]
d.other[[2]]
d.other[[3]]
d.other[[4]]
d.other[[5]]
d.other[[6]]
}
\keyword{Functions}