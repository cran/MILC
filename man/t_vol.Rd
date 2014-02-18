\name{t_vol}
\alias{t_vol}
\title{Predict age at given tumor volume
}
\description{
Function to predict the time (years) at which tumor has reached volume V from the onset of the first malignant cell.           
}
\usage{
t_vol(vol, m)
}
\arguments{
  \item{vol}{tumor volume (mm\eqn{^3})
}
  \item{m}{scale parameter of the Gompertz distribution assumed for the tumor growth
}
}
\value{
Time (years) from the onset of the first malignant cell at which tumor reaches volume 'vol'.
}
\references{
Detterbeck, F. C. and Gibson, C. J. (2008) Turning gray: The natural history of lung cancer over time. \emph{Journal of Thoracic Oncology}, 3(7):781-792.

Laird, A. K. (1964) Dynamics of tumor growth. \emph{British Journal of Cancer}, 18(3):490-502.
}
\author{
Stavroula A. Chrysanthopoulou
}
\seealso{\code{\link{t_prog}}
}
\examples{
# Time (years) elapsed from the onset of the first malignant cell until 
# the tumor reaches volume:

# 0.5 mm^3:
t_vol(0.5, 0.00042)


# 10 mm^3:
t_vol(10, 0.00042)

}
\keyword{Functions}