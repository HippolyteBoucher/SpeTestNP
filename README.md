# SpeTestNP

This R package performs nonparametric tests of parametric specifications. Five heteroskedasticity-robust tests are available: Bierens (1982), Zheng (1996), Escanciano (2006), Lavergne and Patilea (2008), and Lavergne and Patilea (2012). Specific bandwidth and kernel methods can be chosen along with many other options. Allows parallel computing to quickly compute p-values based on the bootstrap. 

Hippolyte Boucher (Hippolyte.Boucher@outlook.com) is the author of `SpeTestNP` and Pascal Lavergne (lavergnetse@gmail.com) is a contributor. Both Hippolyte Boucher and Pascal Lavergne are maintainers and any question or bug should be reported to one of them.

## Installation

To install `SpeTestNP` from CRAN simply run the following command:

    install.packages("SpeTestNP")

To install `SpeTestNP` from Github the package `devtools` should be installed and the following commands should be run:

    library(devtools) 
    install_github("HippolyteBoucher/SpeTestNP")
  
`SpeTestNP` requires the packages `stats`  (already installed and loaded by default in Rstudio), `foreach` and `doParallel` (if parallel computing is used to derive p-values based on the bootstrap) to be installed.

## Examples 

    library(SpeTestNP)
    
    ### Create artificial data, choose a parametric specification, and perform the test
    
    # Example 1: The true model is linear and we estimate the model linearly
    
    n<-500
    x<-rnorm(n)
    y<-x+rnorm(n)*0.7
    
    eq<-lm(y~x+0)
    
    # The test of Bierens (1982) is used with 300 wild bootstrap samples used
    # to compute the p-value
    
    SpeTest(eq=eq,type="icm",rejection="bootstrap",boot="wild",nboot=300)
    
    
    # Example 2: The true model is nonlinear and we estimate the model linearly
    
    n<-500
    x<-rnorm(n)
    y<-x^2+rnorm(n)*0.7
    
    eq<-lm(y~x+0)
    
    # The test of Lavergne and Patilea (2012) is used with a p-value based on the
    # asymptotic distribution of the normalized test statistic under the null hypothesis
    
    summary(SpeTest(eq=eq,type="sicm",rejection="asymptotics"))
    

## References 

H.J. Bierens (1982), ["Conistent Model Specification Test"][1], *Journal of Econometrics*, 20 (1), 105-134

J. C. Escanciano (2006), ["A Consistent Diagnostic Test for Regression Models Using Projections"][3], *Econometric Theory*, 22 (6), 1030-1051

P. Lavergne and V. Patilea (2008), ["Breaking the Curse of Dimensionality in Nonparametric Testing"][4], *Journal of Econometrics*, 143 (1), 103-122

P. Lavergne and V. Patilea (2012), ["One for All and All for One: Regression Checks with Many Regressors"][5], *Journal of Business & Economic Statistics*, 30 (1), 41-52

J.X. Zheng (1996), ["A Consistent Test of Functional Form via Nonparametric Estimation Techniques"][2], *Journal of Econometrics*, 75 (2), 263-289

[1]: https://www.sciencedirect.com/science/article/pii/0304407682901051
[2]: https://econpapers.repec.org/article/eeeeconom/v_3a75_3ay_3a1996_3ai_3a2_3ap_3a263-289.htm
[3]: https://www.jstor.org/stable/4093212
[4]: https://www.sciencedirect.com/science/article/pii/S0304407607001601
[5]: https://www.tandfonline.com/doi/full/10.1198/jbes.2011.07152


