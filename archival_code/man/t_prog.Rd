\name{t_prog}
\alias{t_prog}
\title{Predict age at distinct lung cancer natural history stages
}
\description{
Function to simulate the time points (years) from the onset of the first malignant cell to major states in the natural history of lung cancer, i.e., the time at the beginning of the regional (Treg), and distant (Tdist) stage, and time at diagnosis (Tdiagn). 
}
\usage{
t_prog(N, m, cdiagn, creg, cdist)
}
\arguments{
  \item{N}{Total number of random numbers from the log-Normal distributions
}
  \item{m}{scale parameter of the Gompertz distribution for the tumor growth
}
  \item{cdiagn}{vector with location and scale parameters for the logNormal distribution of the tumor volume at diagnosis.
}
  \item{creg}{vector with location and scale parameters for the logNormal distribution of the tumor volume at the beginning of the regional stage.
}
  \item{cdist}{vector with location and scale parameters for the logNormal distribution of the tumor volume at the beginning of the distant stage.
}
}
\details{
The MILC model employs the Gompertz function to describe the tumor growth, i.e., the tumor volume V(t) at age t is:
\deqn{   \frac{V(t)}{V_0} = e^{\frac{s}{m}(1-e^{-mt})}   }
where, V\eqn{_0} is the minimum tumor volume (one malignant cell), and m, s are the scale and shape parameters of the Gompertz distribution.
Assuming a spherical tumor growth, the tumor size as a function of its diameter d(t) at age t, is:
\deqn{V(t) = \frac {\pi}{6} [d(t)]^3}.

Given a Gompertzian tumor growth, the distribution of tumor volumes at specific time points can be described by a log-Normal distribution, having as starting point the beginning of the local stage (onset of the first malignant cell). We define three different log-Normal distributions to simulate the tumor volume at the beginning of the regional (V\eqn{_{reg}}), and distant stage (V\eqn{_{dist}}), as well as the tumor at diagnosis (V\eqn{_{diagn}}).
}
\value{An R-object of class "list" with the following five components:

[[1]] T_{reg}: time(years) from the onset of the first malignant cell to the beginning of the regional tumor stage

[[2]] T_{dist}: time(years) from the onset of the first malignant cell to the beginning of the distant tumor stage

[[3]] T_{diagn}: time(years) from the onset of the first malignant cell to diagnosis

[[4]] D_{diagn}: tumor diameter at diagnosis

[[5]] stage: tumor stage at diagnosis
}
\references{
Detterbeck, F. C. and Gibson, C. J. (2008) Turning gray: The natural history of lung cancer over time. \emph{Journal of Thoracic Oncology}, 3(7):781-792.

Laird, A. K. (1964) Dynamics of tumor growth. \emph{British Journal of Cancer}, 18(3):490-502.

Meza, R., Hazelton, W. D., Colditz, G. A. and Moolgavkar, S. H. (2008) Analysis of lung cancer incidence in the nurses' health and the health professionals' follow-up studies using a multistage carcinogenesis model. \emph{Cancer Causes & Control}, 19(3):317-328.

McMahon, P. M. (2005) Policy assessment of medical imaging utilization: methods and applications [doctoral thesis]. PhD thesis.

Koscielny, S., Tubiana, M. and Valleron, A. J. (1985) A simulation model of the natural history of human breast cancer. \emph{Br J Cancer}, 52(4):515-524.

}
\author{
Stavroula A. Chrysanthopoulou
}
\seealso{\code{\link{t_vol}}
}
\examples{
# Predicted disease progression given certain values of the log-normal and Gompertz 
# distributions assumed for the tumor volumes and growth in the MILC model.
 t_prog (1, 0.00042, c(3.91, 3.91), c(1.1, 1.1), c(2.8, 2.8))
}
\keyword{Functions}