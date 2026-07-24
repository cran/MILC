\name{HT_mal}
\alias{HT_mal}
\title{Calculate cumulative hazard for the 1st malignant cell
}
\description{
Cumulative hazard up to age t for the onset of first malignant cell, given gender, smoking status and smoking intensity.
}
\usage{
HT_mal(t, gender, d, smok_status)
}
\arguments{
  \item{t}{ age
}
  \item{gender}{ gender ("male", or "female")
}
  \item{d}{ smoking instensity (number of cigarettes)
}
  \item{smok_status}{ whether the individual is a smoker ("Y") or not ("N")
}
}
\details{Integration from age 0 to age t of the hazard rate h(t) for the onset of the first malignant cell, as described by the Two-Stage Clonal Expansion (TSCE) carcinogenesis model:

\deqn{
h(t) = \frac{\nu \mu X (e^{(\gamma+2B)t}-1)}{\gamma+B(e^{(\gamma+2B)t}+1)}
     }

with \eqn{\gamma = \alpha - \beta -\mu} and B=\eqn{\frac{1}{2}(-\gamma + \sqrt{\gamma^2 +4 \alpha \mu})}

Power laws express the effect of smoking on the hazard for the onset of the first malignant cell, given smoking intensity d at age t, i.e.:

\tabular{rcl}{
\eqn{\alpha = \alpha_0 [1+\alpha_1 \cdot d^{\alpha_2}]} \tab and \tab \eqn{\gamma = \gamma_0 [1+\alpha_1 \cdot d^{\alpha_2}] }
}

where

\tabular{rl}{
X:	     	\tab   total number of normal stem cells	\cr
\eqn{\nu}:    \tab   normal cell initiation rate		\cr
\eqn{\alpha}: \tab   division rate of initiated cells 		\cr
\eqn{\beta}:  \tab   apoptosis rate of initiated cells		\cr
\eqn{\mu}:    \tab   malignant conversion rate of initiated cells \cr
d:	      \tab  average number of cigarettes smoked per day
	}

}
\value{ Cumulative hazard for the onset of the first malignant cell.
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
\seealso{ \code{\link{ht_mal_int}, \link{t_mal}}
}
\examples{
# Cumulative hazard of the onset of the first malingant cell for a man, 50 years old, 
# who smokes 15 cigarettes per day on average.
HT_mal(50, "male", 15, "Y")
}
\keyword{Functions}