\name{SpeTestNP}
\alias{SpeTestNP}

\keyword{test}
\keyword{nonparametric}
\keyword{parametric}

\concept{nonparametric test}
\concept{specification test}
\concept{linear model}
\concept{parametric model}


\title{
Nonparametric specification test package
}

\description{
This package performs nonparametric tests of parametric specifications. Five heteroskedasticity-robust tests are available: Bierens (1982), Zheng (1996), Escanciano (2006), Lavergne and Patilea (2008), and Lavergne and Patilea (2012). Specific bandwidth and kernel methods can be chosen along with many other options. Allows parallel computing to quickly compute p-values based on the bootstrap. 
}

\details{
The package contains the following functions:
  \itemize{
  \item \code{\link{SpeTest}} returns a specfication test object of class \code{STNP}
  \item \code{\link{print.STNP}} a method to print a specification test \code{STNP} object
  \item \code{\link{summary.STNP}} a method to print a summary of a specification test \code{STNP} object
  \item \code{\link{SpeTest_Stat}} computes the test statistic
  \item \code{\link{SpeTest_Dist}} generates a vector drawn from the distribution of the test statistic under the null hypothesis using the bootstrap
  }
}

\note{
  Requires the packages \code{stats}  (already installed and loaded by default in Rstudio), \code{foreach}, \code{parallel} and \code{doParallel} (if parallel computing is used to generate the vector) to be installed

  To use the package to its full potential see the references
}

\author{
Hippolyte Boucher <Hippolyte.Boucher@outlook.com>

Pascal Lavergne <lavergnetse@gmail.com>
}

\references{
H.J. Bierens (1982), "Consistent Model Specification Test",
\cite{Journal of Econometrics}, 20 (1), 105-134

J.C. Escanciano (2006), "A Consistent Diagnostic Test for Regression Models using Projections",
\cite{Economic Theory}, 22 (6), 1030-1051

P. Lavergne and V. Patilea (2008), "Breaking the Curse of Dimensionality in Nonparametric Testing", \cite{Journal of Econometrics}, 143 (1), 103-122

P. Lavergne and V. Patilea (2012), "One for All and All for One: Regression Checks with Many Regressors",
\cite{Journal of Business and Economic Statistics}, 30 (1), 41-52

J.X. Zheng (1996), "A Consistent Test of Functional Form via Nonparametric Estimation Techniques",
\cite{Journal of Econometrics}, 75 (2), 263-289
}
