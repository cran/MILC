\name{nat_hist}
\alias{nat_hist}
\title{Simulate one individual trajectory
}
\description{Function to run the MILC model and predict a full lung cancer trajectory depending on the age, gender and smoking history of an individual (see \code{\link{MILC-package}} for more details).
}
\usage{
nat_hist(dat, pred_yrs, gender, status, ts, tq, m, cdiagn, creg, cdist)
}
\arguments{
  \item{dat}{7-dimensional numeric vector, with the first 5 values being random numbers from Unif[0,1], required for the simulation, and the 6th and 7th value referring to the age (years) and smoking intensity (number of cigarettes) respectively.
}
  \item{pred_yrs}{prediction period (years)
}
  \item{gender}{ "male" or "female"
}
  \item{status}{smoking status, i.e., whether the person is "never", "former", or "current" smoker
}
  \item{ts}{start smoking age (years), when relevant
}
  \item{tq}{quit smoking age (years), when relevant
}
  \item{m}{scale parameter of the Gompertz distribution assumed for the tumor growth
}
  \item{cdiagn}{2-dimensional vector with location and scale parameters of the log-Normal distribution assumed for the tumor volume at diagnosis
}
  \item{creg}{2-dimensional vector with location and scale parameters of the log-Normal distribution assumed for the tumor volume at the beginning of the regional stage

}
  \item{cdist}{2-dimensional vector with location and scale parameters of the log-Normal distribution assumed for the tumor volume at the beginning of the distant stage	}
}
\value{An R-object of type "list" with the following 20 components:

[[1]]: T-entry: age(years) at the beginning of the prediction period

[[2]]: T_mal: age(years) at the onset of the first malignant cell

[[3]]: T_reg: age(years) at the beginning of the regional stage

[[4]]: T_dist: age(years) at the beginning of the distant stage

[[5]]: T_diagn: age(years) at diagnosis

[[6]]: D_diagn: tumor diameter at diagnosis

[[7]]: stage: tumor stage at diagnosis

[[8]]: T-pred: age(years) at the end of the prediction period

[[9]]: T_do: predicted age(years) at death from a cause other than lung cancer

[[10]]: T_dl: predicted age(years) at death from lung cancer

[[11]]: T_final: age (years) at the end of the simulated trajectory

[[12]]: lung_inc: whether the person developed (1="Yes") lung cancer or not (0="No")

[[13]]: excl: exclude unreasonable cases ("Yes", "No")

[[14]]: cause: cause of death ("lung", "other", NA)

[[15]]: T_death: age(years) at death from any cause

[[16]]: gender

[[17]]: smoking status

[[18]]: start smoking age(years) 

[[19]]: quit smoking age(years) 

[[20]]: smoking intensity (number of cigarettes)
}
\author{
Stavroula A. Chrysanthopoulou}
\seealso{\code{\link{t_mal}, \link{t_prog}, \link{tdeath_other}, \link{tdeath_lung}}
}
\examples{
# In the following examples we predict lung cancer trajectories for a man, 50 years old 
# at the beginning of the prediction period, who has started smoking at the age of 20 years
# and smokes 30 cigerettes per day on average. The model predicts 20 years ahead. 

# We present three possible trajectories:

# In the first case the person does not die before the end of the prediction period:
set.seed(33)
nat_hist ( c(runif(5),50,30), 20, "male", "current", 20, NA, 
		0.00042, c(3.91, 3.91), c(1.1, 1.1), c(2.8, 2.8))

# In the second case the person dies at the age of 62.43 years from lung cancer:
set.seed(1470)
nat_hist ( c(runif(5),50,30), 20, "male", "current", 20, NA, 
		0.00042, c(3.91, 3.91), c(1.1, 1.1), c(2.8, 2.8))

# In the third case the person dies at the age of 69.68 years from a cause
# other than lung cancer:
set.seed(1450)
nat_hist ( c(runif(5),50,30), 20, "male", "current", 20, NA, 
		0.00042, c(3.91, 3.91), c(1.1, 1.1), c(2.8, 2.8))
}
\keyword{Functions}