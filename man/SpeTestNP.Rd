\name{SpeTestNP}
\alias{SpeTestNP}

\keyword{test}
\keyword{nonparametric}

\concept{nonparametric test}
\concept{specification test}
\concept{linear model}


\title{
Nonparametric Specification Tests for Linear Models
}

\description{
\pkg{SpeTestNP} implements various heteroskedasticity robust nonparametric tests of specification for linear models and non-linear additive models: Bierens (1982) Integrated Conditional Moment test; Zheng (1996) test; Escanciano (2006) test; Lavergne and Patilea (2008); Lavergne and Patilea (2012) Smooth Integrated Conditional Moment test.
}

}
\details{
The package contains the following functions:
  \itemize{
  \item \code{\link{SpeTest}} returns an object of class STNP
  \item \code{\link{SpeTest_Stat}} computes the test statistic of choice
  \item \code{\link{SpeTest_Dist}} computes the distribution of the test statistic of choice
  }
}

\note{
  Requires the packages "np" (version >= 0.60), "pracma"

  To use the package to its full potential see the references.
}

\author{
Hippolyte Boucher <Hippolyte.Bocuher@outlook.com>

Maintainer: Pascal Lavergne <lavergnetse@gmail.com>
}

\references{
H.J. Bierens (1982) Consistent Model Specification Test,
\cite{Journal of Econometrics}, 20, 105-134

J.C. Escanciano (2006) A Consistent Diagnostic Test for Regression Models using Projections,
\cite{Economic Theory}, 22, 1030-1051

Gozalo, P.L., (1997) Nonparametric bootstrap analysis with applications to demographic effects in demand functions,
\cite{Journal of Econometrics}, 81, 357-393.

P. Lavergne and V. Patilea (2008) Breaking the urse of dimensionality in nonparametric testing,
\cite{Journal of Econometrics}, 143, 103-

P. Lavergne and V. Patilea (2012) One for All and All for One: Regression Checks with Many Regressors,
\cite{Journal of Business and Economic Statistics}, 30, 41-52

Li, Q., S. Wang (1998) A simple consistent bootstrap test for a parametric
regression function,
\cite{Journal of Econometrics}, 87, 145-165

Yin, J., Z. Geng, R. Li, H. Wang (2010) Nonparametric covariance model,
\cite{Statistica Sinica}, 20, 469-479

J.X. Zheng (1996) A consistent test of functional form via nonparametric estimation techniques,
\cite{Journal of Econometrics}, 75, 263-289

}
