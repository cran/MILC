\name{ht_mal_int}
\alias{ht_mal_int}
\title{Integrating hazard for the onset of th 1st malignant cell
}
\description{Cumulative hazard for the onset of the first malignant cell between two time points (ages).
}
\usage{
ht_mal_int(lower, upper, g, d, smoking)
}
\arguments{
  \item{lower}{ age (years), lower limit of the integral
}
  \item{upper}{ age (years), upper limit of the integral
}
  \item{g}{ gender ("male", or "female")
}
  \item{d}{ smoking intensity group based on the average number of cigarettes smoked per day
}
  \item{smoking}{ binary variable indicating whether the person was smoking between ages "lower" and "upper" ("Y"=yes, "N"=no)
}
}
\details{Integration from age \emph{`lower'} to age \emph{`upper'} of the hazard rate h(t) for the onset of the first malignant cell, as described by the Two-Stage Clonal Expansion (TSCE) carcinogenesis model (see function '\code{\link{HT_mal}}' for more details).
}
\value{The value of the calculated integral, equal to the probability of developing the first malignant cell between ages "lower" and "upper".
}
\references{
Moolgavkar, S. H.,  and Luebeck, G. (1990) Two-event model for carcinogenesis: Biological, mathematical, and statistical considerations. \emph{Risk Analysis}, 10(2):323-341.

Heidenreich, W. F., Luebeck, E. G. and Moolgavkar, S. H. (1997) Some properties of the hazard function of the two-mutation clonal expansion model. \emph{Risk Analysis}, 17(3):391-399.

Hazelton, W. D., Clements, M. S. and Moolgavkar, S. H. (2005) Multistage carcinogenesis and lung cancer mortality in three cohorts. \emph{Cancer Epidemiology Biomarkers & Prevention}, 14(5):1171-1181.

Hazelton, W. D., Luebeck,  E. G., Heidenreich, W. E. and Moolgavkar, S. H.(2001) Analysis of a historical cohort of chinese tin miners with arsenic, radon, cigarette smoke, and pipe smoke exposures using the biologically based two-stage clonal expansion model. \emph{Radiation Research}, 156(1):78-94.

Meza, R., Hazelton, W. D., Colditz, G. A. and Moolgavkar, S. H. (2008) Analysis of lung cancer incidence in the nurses' health and the health professionals' follow-up studies using a multistage carcinogenesis model. \emph{Cancer Causes & Control}, 19(3):317-328.
}
\author{
Stavroula A. Chrysanthopoulou}
\seealso{\code{\link{HT_mal}, \link{t_mal}}
}
\examples{
# The following examples calculate the probability of a man, who on average smokes 
# 30 cigarettes per day, to develop the first malignant cell:
ht_mal_int(0, 50, "male", 30, "Y")	# before the age of 50
ht_mal_int(0, 75, "male", 30, "Y")	# before the age of 75
ht_mal_int(45, 75, "male", 30, "Y")	# between the ages of 45 and 75
}
\keyword{Functions}