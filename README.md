# SpeTestNP

This R package performs nonparametric test of linearity. Five heteroskedasticity-robust tests are available: Bierens (1982), Zheng (1996), Escanciano (2006), Lavergne and Patilea (2008), and Lavergne and Patilea (2012). Specific bandwidth and kernel methods can be chosen along with many other options.

To install the package from Github the package `devtools` should be installed and the following commands should be run:

    library(devtools) 
    install_github("HippolyteBoucher/SpeTestNP")
  


## Examples 

    library(SpeTestNP)
    
    ### Create artificial data and perform the test
    
    n<-500
    x<-rnorm(n)
    
    # The true model is linear and the test of Zheng (1996) is used with 300 wild bootstrap samples
    
    y<-x+rnorm(n)*0.7
    
    eq<-lm(y~x)
    SpeTest(eq,type="zheng",boot="wild",nboot=300)
    
    # The true model is nonlinear and the test of Lavergne and Patilea (2012) is used
    # with 99 smooth conditional moments bootstrap samples, with 100 values in the
    # hypersphere and starting values equal to 0, 1, and -1, and the result is summarized
    
    y<-x^2+rnorm(n)*0.7
    summary(SpeTest(x=x,y=y,inter=T,type="sicm",boot="smooth",
    nbeta=100,direct=c(0,1,-1)))
    



## References 

H.J. Bierens, ["Conistent Model Specification Test"][1], *Journal of Econometrics*, 20 (1), 105-134

J.X. Zheng, ["A Consistent Test of Functional Form via Nonparametric Estimation Techniques"][2], *Journal of Econometrics*, 75 (2), 263-289

J. C. Escanciano, ["A Consistent Diagnostic Test for Regression Models Using Projections"][3], *Econometric Theory*, 22 (6), 1030-1051

P. Lavergne and V. Patilea, ["Breaking the Curse of Dimensionality in Nonparametric Testing"][4], *Journal of Econometrics*, 143 (1), 103-122

P. Lavergne and V. Patilea, ["One for All and All for One: Regression Checks with Many Regressors"][5], *Journal of Business & Economic Statistics*, 30 (1), 41-52

[1]: https://www.sciencedirect.com/science/article/pii/0304407682901051
[2]: https://econpapers.repec.org/article/eeeeconom/v_3a75_3ay_3a1996_3ai_3a2_3ap_3a263-289.htm
[3]: https://www.jstor.org/stable/4093212
[4]: https://www.sciencedirect.com/science/article/pii/S0304407607001601
[5]: https://www.tandfonline.com/doi/full/10.1198/jbes.2011.07152


