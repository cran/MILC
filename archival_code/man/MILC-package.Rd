\name{MILC-package}
\alias{MILC-package}
\alias{MILC}
\docType{package}
\title{
MIcrosimulation Lung Cancer (MILC) model
}
\description{
The MILC package is designed to predict individual trajectories describing the onset of, progression, and (potentially) death from lung cancer using the MILC continuous time microsimulation model. 
}
\details{
The MIcrosimulation Lung Cancer model (MILC), is a new, dynamic, continuous time microsimulation model that, in its current version, comprises a module that describes the natural history of lung cancer in the absence of any screening or treatment component. The model simulates the course of lung cancer from the disease free state to the local, regional, and distant states and eventually to death from either lung cancer or some other cause. When predicting individual trajectories, the model accounts for age, gender, and smoking history, including smoking status, start and quit smoking ages, and average number of cigarettes smoked per day when relevant.

The model comprises four main components:
\itemize{
\item \bold{Onset of the first malignant cell:} The local stage of the lung cancer tumor initiates with the onset of the first malignant cell, as described by the Two-Stage Clonal Expansion (TSCE) carcinogenesis model (see \code{\link{HT_mal}} for more details).
\item \bold{Tumor growth:} The model assumes a spherical tumor growth described by a Gompertz distribution (see \code{\link{t_prog}} for more details).
\item \bold{Disease progression:} Given a Gompertzian tumor growth, the tumor volume at specific time points is described by log-Normal distributions (see \code{\link{t_prog}} for more details).
\item \bold{Survival:} The model employs the Cumulative Incidence Function (CIF) non-parametric technique to simulate survival in a competing-risks frame.
}

\tabular{ll}{
Package: \tab MILC\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-02-18\cr
License: \tab GPL-2\cr
}
}
\author{
Stavroula A. Chrysanthopoulou <Stavroula_Chrysanthopoulou@brown.edu>
}
\references{
Rutter, C. M., Zaslavsky, A. M. and Feuer, E.J. (2011) Dynamic microsimulation models for health outcomes. \emph{Medical Decision Making}, 31(1):10-18.

Moolgavkar, S. H.,  and Luebeck, G. (1990) Two-event model for carcinogenesis: Biological, mathematical, and statistical considerations. \emph{Risk Analysis}, 10(2):323-341.

Heidenreich, W. F., Luebeck, E. G. and Moolgavkar, S. H. (1997) Some properties of the hazard function of the two-mutation clonal expansion model. \emph{Risk Analysis}, 17(3):391-399.

Hazelton, W. D., Clements, M. S. and Moolgavkar, S. H. (2005) Multistage carcinogenesis and lung cancer mortality in three cohorts. \emph{Cancer Epidemiology Biomarkers & Prevention}, 14(5):1171-1181.

Hazelton, W. D., Luebeck,  E. G., Heidenreich, W. E. and Moolgavkar, S. H.(2001) Analysis of a historical cohort of chinese tin miners with arsenic, radon, cigarette smoke, and pipe smoke exposures using the biologically based two-stage clonal expansion model. \emph{Radiation Research}, 156(1):78-94.

Meza, R., Hazelton, W. D., Colditz, G. A. and Moolgavkar, S. H. (2008) Analysis of lung cancer incidence in the nurses' health and the health professionals' follow-up studies using a multistage carcinogenesis model. \emph{Cancer Causes & Control}, 19(3):317-328.

Detterbeck, F. C. and Gibson, C. J. (2008) Turning gray: The natural history of lung cancer over time. \emph{Journal of Thoracic Oncology}, 3(7):781-792.

Laird, A. K. (1964) Dynamics of tumor growth. \emph{British Journal of Cancer}, 18(3):490-502.

McMahon, P. M. (2005) Policy assessment of medical imaging utilization: methods and applications [doctoral thesis]. PhD thesis.

Koscielny, S., Tubiana, M. and Valleron, A. J. (1985) A simulation model of the natural history of human breast cancer. \emph{Br J Cancer}, 52(4):515-524.

Klein, J. P., and Moeschberger, M. L. (2003) Survival analysis: techniques for censored and truncated data.

Gray, B. (2013) cmprsk: Subdistribution Analysis of Competing Risks, URL: http://CRAN.R-project.org/package=cmprsk. R package version 2.2-6.
}

