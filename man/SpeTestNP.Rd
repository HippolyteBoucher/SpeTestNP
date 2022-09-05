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
\pkg{SpeTestNP} performs nonparametric test of linearity. Five heteroskedasticity-robust tests are available: Bierens (1982), Zheng (1996), Escanciano (2006), Lavergne and Patilea (2008), and Lavergne and Patilea (2012). Specific bandwidth and kernel methods can be chosen along with many other options.
}

}
\details{
The package contains the following functions:
  \itemize{
  \item \code{\link{SpeTest}} returns an object of class STNP
  \item \code{\link{SpeTest_Stat}} computes the test statistic of choice
  \item \code{\link{SpeTest_Dist}} computes the distribution of the test statistic of choice under the null
  }
}

\note{
  Requires the packages "np" (version >= 0.60), "pracma"

  To use the package to its full potential see the references
}

\author{
Hippolyte Boucher <Hippolyte.Bocuher@outlook.com>

Maintainer: Pascal Lavergne <lavergnetse@gmail.com>
}

\references{
H.J. Bierens (1982), "Consistent Model Specification Test",
\cite{Journal of Econometrics}, 20 (1), 105-134

J.C. Escanciano (2006), "A Consistent Diagnostic Test for Regression Models using Projections",
\cite{Economic Theory}, 22 (6), 1030-1051

Gozalo, P.L., (1997), "Nonparametric Bootstrap Analysis with Applications to Demographic Effects in Demand Functions",
\cite{Journal of Econometrics}, 81 (2), 357-393.

P. Lavergne and V. Patilea (2008), "Breaking the Curse of Dimensionality in Nonparametric Testing",
\cite{Journal of Econometrics}, 143 (1), 103-122

P. Lavergne and V. Patilea (2012), "One for All and All for One: Regression Checks with Many Regressors",
\cite{Journal of Business and Economic Statistics}, 30 (1), 41-52

Li, Q., S. Wang (1998), "A Simple Consistent Bootstrap Test for a Parametric
Regression Function",
\cite{Journal of Econometrics}, 87 (1), 145-165

Yin, J., Z. Geng, R. Li, H. Wang (2010), "Nonparametric covariance model",
\cite{Statistica Sinica}, 20 (1), 469-479

J.X. Zheng (1996), "A Consistent Test of Functional Form via Nonparametric Estimation Techniques",
\cite{Journal of Econometrics}, 75 (2), 263-289

}
