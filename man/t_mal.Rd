\name{t_mal}
\alias{t_mal}
\title{Predict age at the onset of the first malignant cell
}
\description{Function to predict the age at the onset of the first malignant given gender and smoking history, when relevant.
}
\usage{
t_mal(u, g, ts, tq, d)
}
\arguments{
  \item{u}{random number required for the prediction}
  \item{g}{gender ("male" or "female")}
  \item{ts}{start smoking age (years)}
  \item{tq}{quit smoking age (years)}
  \item{d}{smoking intensity expressed as average number of cigarettes smoked per day}
}
\details{The age at the onset for the first malignant cell is simulated based on the relevant hazard h(t) described by the TSCE carcinogenesis model (see function '\code{\link{HT_mal}}' for more details). Smoking history includes start and quit smoking ages, and smoking intensity expressed as average number of cigarettes smoked per day.
}
\value{Predicted age (years) at which a person may develop the first malignant cell.
}
\references{
Moolgavkar, S. H.,  and Luebeck, G. (1990) Two-event model for carcinogenesis: Biological, mathematical, and statistical considerations. \emph{Risk Analysis}, 10(2):323-341.

Heidenreich, W. F., Luebeck, E. G. and Moolgavkar, S. H. (1997) Some properties of the hazard function of the two-mutation clonal expansion model. \emph{Risk Analysis}, 17(3):391-399.

Hazelton, W. D., Clements, M. S. and Moolgavkar, S. H. (2005) Multistage carcinogenesis and lung cancer mortality in three cohorts. \emph{Cancer Epidemiology Biomarkers & Prevention}, 14(5):1171-1181.

Hazelton, W. D., Luebeck,  E. G., Heidenreich, W. E. and Moolgavkar, S. H.(2001) Analysis of a historical cohort of chinese tin miners with arsenic, radon, cigarette smoke, and pipe smoke exposures using the biologically based two-stage clonal expansion model. \emph{Radiation Research}, 156(1):78-94.

Meza, R., Hazelton, W. D., Colditz, G. A. and Moolgavkar, S. H. (2008) Analysis of lung cancer incidence in the nurses' health and the health professionals' follow-up studies using a multistage carcinogenesis model. \emph{Cancer Causes & Control}, 19(3):317-328.
}
\author{Stavroula A. Chrysanthopoulou}
\seealso{ \code{\link{HT_mal}, \link{ht_mal_int}, \link{nat_hist}}
}
\examples{
# The following examples predict the possible ages at the onset of the first 
# malignant cell of:

set.seed(7461) 	
# a man who has never smoked:
   t_mal(runif(1), "male", NA, NA, NA) 

# a man, current smoker,who started smoking at 20 years:
   t_mal(runif(1), "male", 20, NA, 25)

# a man, former smoker, who started and quit smoking at 20 and 50 years respectively, 
# and used to smoke 25 cigarettes per day on average:
   t_mal(runif(1), "male", 20, 50, 25) 
}
\keyword{Functions}