\name{SpeTest}
\alias{SpeTest}

\keyword{test}
\keyword{nonparametric}

\concept{nonparametric test}
\concept{specification test}
\concept{linear model}

\title{
Nonparametric Specification Tests
}

\description{
\code{SpeTest} returns the statistic and the bootstrap p-value for five different non-parametric specification tests of linearity of the model. The linear specification of the model is tested against all alternatives which are not too small, typically these alternatives are bigger than the square-root of the sample size but it depends on the test which is used.
}


\usage{
SpeTest(eq, x, y, inter=TRUE, type="icm", boot="wild", nboot=99, ker="normal",
knorm="sd", hv="default", cch="default", nbeta=200, beta0="default", 
direct="default")
}


\arguments{
  \item{eq}{
An object of \code{\link{class}} \code{"\link{lm}"} (uneeded if x and y are provided)
}
\item{x}{
A matrix (whose columns are not perfect linear combinations of each other) or a vector of regressors (uneeded if eq is provided)
}
\item{y}{
A vector of dependent variables (uneeded if eq is provided)
}
\item{inter}{
In case if x and y are provided an intercept is included if \code{inter = TRUE}
}
\item{type}{
To use Bierens (1982) Integrated Conditional Moment test \code{type = "icm"}; To use Zheng (1996) test \code{type = "zheng"}; To use Escanciano (2006) test \code{type = "esca"}; To use Lavergne and Patilea (2008) test \code{type = "pala"}; To use Lavergne and Patilea (2012) Smooth Integrated Conditional Moment test \code{type = "sicm"}
}
\item{boot}{
Bootstraps method, \code{boot = "wild"} to use wild bootstraps (Li, Wang 1998),
\code{boot = "smooth"} to use smooth conditional moments bootstraps (Gozalo 1997)
}
\item{nboot}{
Number of bootstraps used to build the p-value, by default \code{nboot = 99}
}
\item{ker}{
Kernel function used in the central matrix,

If \code{ker = "normal"} the kernel function is the normal p.d.f;
If \code{ker = "triangle"} the kernel function is the triangular p.d.f;
If \code{ker = "logistic"} the kernel function is the logistic p.d.f;
If \code{ker = "sinc"} the kernel function is the sine cardinal function.
}
\item{knorm}{
Normalization of the kernel,

If \code{knorm = "sd"} then the standard deviation using the kernel function equals 1; If \code{knorm ="sq"} then the integral of the squared kernel function equals 1.
}
\item{hv}{
Non-parametric errors covariance estimator bandwidth

By \code{"default"} the bandwidth is the scaled rule of thumb \code{cch = 1.06*n^(-1/(4+k))}.

The user may also give a specific value for the bandwidth.
}
\item{cch}{
Central kernel matrix bandwidth,

If \code{type = "icm"} or \code{type = "esca"}, \code{cch} always equals \code{1}, hence changing the bandwidth won't affect the statistic.

For \code{type = zheng} the \code{"default"} bandwidth is the scaled
rule of thumb: \code{cch = 1.06*n^(-1/5)}.

For \code{type = "sicm"} and \code{type = "pala"} the \code{"default"}
bandwidth is the scaled rule of thumb: \code{cch = 1.06*n^(-1/(4+k))}

The user may change the bandwidth for \code{"zheng"}, \code{"sicm"} and \code{"pala"}.
}
\item{nbeta}{
If the Patilea and Lavergne or the Smooth Integrated Conditional Moment test
are used, \code{nbeta} is the number of "betas" in the unit hypersphere used to compute the statistic. Computing time increases with \code{nbeta} large.

By \code{"default"} it is the \code{200} times the square root of the regressors dimension.

Do not modify unless you have extensive knowledge about the test.
}
\item{beta0}{
If the Patilea and Lavergne test is used, \code{beta0} is the initial beta from which
the statistic is computed, by \code{"default"} it is the OLS estimator.

Do not modify unless you have extensive knowledge about the test.
}

\item{direct}{
If the Smooth Integrated Conditional Moment test is used, \code{direct} is the
initial prediction for beta; This direction should be a vector of \code{0} (for
no direction), \code{1} (for positive direction) and \code{-1} (for negative direction)

By \code{"default"} no direction are given to the hypersphere.
Do not modify unless you have extensive knowledge about the test.

For ex, \code{c(1,-1,0)} indicates that the user thinks that the 1st regressor
has a positive effect on the dependent variable, that the 2nd regressor
has a negative effect on the dependent variable, and that he has no idea
about the effect of the 3rd regressor.

Note that giving a direction to the intercept does not modify the statistic therefore
the 1st direction in your vector is given to the 1st regressor which is not
the intercept.
}
}


\details{
By default to perform the test of Bierens (1982), only an object of class \code{"lm"} or a dependent variable \code{y} and a set of covariates \code{x} are needed.

Note that if both an \code{"lm"} model and dependent variable \code{y} and regressors \code{x} are provided or if none are provided then the test will not be performed.

Other tests are available: Zheng (1996), Escanciano (2006), Lavergne and Patilea (2008) and Lavergne and Patilea Smooth Integrated Conditional Moments test (2012).

A lot more things can and should be specified: the bandwidths; the bootstraps' type and the number of bootstraps samples; the kernel functions and the standardization; if the user has knowledge of the tests coined by Lavergne and Patilea he may choose a higher number of betas for the hypersphere (which may slow significantely the process) and an initial "direction" to the hypersphere for the SICM test (none is given by \code{"default"}) or a starting beta for the PALA test (which is the OLS estimator by \code{"default"}).

The nonparametric conditional variance estimator used in the tests is discussed in Yin, Geng, Li, Wang (2010). The wild bootstraps of Li and Wang (1997) and the smooth conditional moments bootstraps of Gozalo (1997) can be used to compute the test p-values.
}


\value{
\code{SpeTest} returns an object of \code{\link{class}} \code{"STNP"}.

\code{\link{summary}} and \code{\link{print}} can be used on objects of this class.

An object of class \code{"STNP"} is a list which contains at least the following elements:

\item{stat}{The value of the statistic used in the test}
\item{pval}{The bootstrap p-value}
\item{type}{The test's type}
\item{nboot}{The number of bootstrap samples used to compute the p-value}
\item{ker}{The kernel function: \code{"normal"} for the normal p.d.f, \code{"triangle"}
for the triangular p.d.f, \code{"logistic"} for the logistic p.d.f and
\code{"sinc"} for the sine cardinal function}
\item{knorm}{The kernel matrix standardisation: \code{"sq"} if the second moment equals
1 or \code{"sd"} if the standard deviation equals 1
}
\item{boot}{The bootstraps' type}
}


\note{
\code{summary.STNP} gives a summary of the test with all the options used.

Requires the packages "foreach", "pracma" and "np".

For more information and to be able to use the package to its full potential see the references.
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


\examples{
\code{
n <- 100
k <- 2
x <- matrix(rnorm(n*k),ncol=k)
y<-1+x\%*\%(1:k)+rnorm(n)
eq<-lm(y~x+0)

SpeTest(x=x,y=y,type="zheng",boot="wild",nboot=300)

summary(SpeTest(x=x,y=y,inter=T,type="sicm",boot="smooth",
nbeta=100,direct=c(0,1,-1)))
}
}
