# APPEstimation
R function to calculate model performance measure adjusted for predictor distributions.

This package provides the function to estimate model performance measures (L1, L2, C-statistics). The difference in the distribution of predictors between two datasets (training and validation) is
adjusted by a density ratio estimate.

"appe.glm" calculates adjusted C statistics by predictor distributions for a generalized linear model with binary
outcome.

"appe.lm" calculates adjusted L1 and L2 errors by predictor distributions for a linear model.


Examples:

n0 = 100

Z = cbind(rnorm(n0,0), rnorm(n0,0))

Y = apply(Z, 1, function(xx) { rlnorm(1, sum(c(0,1,1) * c(1,xx))) })

dat = data.frame(Za=Z[,1], Zb=Z[,2], Y=Y)

mdl = lm(Y~ Za + Zb, data=dat)

n1 = 100

Z1 = cbind(rnorm(n1,-0.5), rnorm(n1,0.5))

Y1 = apply(Z1, 1, function(xx) { rlnorm(1, sum(c(0,1,1) * c(1,xx))) })

dat1 = data.frame(Za=Z1[,1], Zb=Z1[,2], Y=Y1)

appe.lm(mdl, dat, dat1, reps=0)
