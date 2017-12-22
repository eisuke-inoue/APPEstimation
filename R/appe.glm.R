#' C-statistics adjusted for predictor distributions
#' Calculates adjusted \eqn{C} statistics by predictor distributions for a generalized linear model with binary outcome.
#' @export
#' @import stats
#' @param mdl A glm object describing a prediction model to be evaluated.
#' @param dat.train A dataframe used to construct a prediction model (specified in mdl), corresponding to a training data. Need to include outcome and all predictors.
#' @param dat.test A dataframe corresponding to a validation (testing) data. Need to include outcome and all predictors.
#' @param method uLSIF or KLIEP. Same as the argument in densratio function from densratio package.
#' @param sigma A positive numeric vector corresponding to candidate values of a bandwidth for Gaussian kernel. Same as the argument in densratio function from densratio package.
#' @param lambda A positive numeric vector corresponding to candidate values of a regularization parameter. Same as the argument in densratio function from densratio package.
#' @param kernel_num A positive integer corresponding to number of kernels. Same as the argument in densratio function from densratio package.
#' @param fold A positive integer corresponding to a number of the folds of cross-validation in the KLIEP method. Same as the argument in densratio function from densratio package.
#' @param stabilize A logical value as to whether tail weight stabilization is performed or not. If TRUE, both tails of the estimated density ratio distribution are replaced by the constant value which is specified at qstb option.
#' @param qstb A positive numerical value less than 1 to control the degree of weight stabilization. Default value is 0.025, indicating estimated density ratio values less than the 2.5 percentile and more than the 97.5 percentile are set to 2.5 percentile and 97.5 percentile, respectively.
#' @param reps A positive integer to specify bootstrap repetitions. If 0, bootstrap calculations are not performed.
#' @param conf.level A numerical value indicating a confidence level of interval.
#' @return matrix
appe.glm <-
function(mdl, dat.train, dat.test, method="uLSIF", sigma=NULL, lambda=NULL, kernel_num=NULL, fold=5, stabilize=TRUE, qstb=0.025, reps=2000, conf.level=0.95) {

    n0 = nrow(dat.train)
    n1 = nrow(dat.test)
    on = as.character(formula(mdl$call)[[2]])

    ## observed & predicted response values
    Y1   = dat.test[,on]
    scr1 = predict(mdl, newdata=dat.test)
    scr0 = predict(mdl, newdata=dat.train)

    ## weight calculation via package 'densratio'
    xtrain = update(mdl, data=dat.train, x=TRUE)$x[,-1,drop=FALSE]
    xtest  = update(mdl, data=dat.test,  x=TRUE)$x[,-1,drop=FALSE]
    wgt1   = densratio.appe(scr0,   scr1,  method, sigma, lambda, kernel_num, fold, stabilize, qstb)
    wgt2   = densratio.appe(xtrain, xtest, method, sigma, lambda, kernel_num, fold, stabilize, qstb)

    ## predictive performance measure
    Cv   = cvalest.bin(Y1, scr1)
    Cvw1 = cvalest.bin(Y1, scr1, wgt1)
    Cvw2 = cvalest.bin(Y1, scr1, wgt2)

    message("\nPoint estimates:")
    result = data.frame(c(Cv, Cvw1, Cvw2))
    names(result) = 'Est'
    row.names(result) = c('Cstat','C adjusted by score','C adjusted by predictors')
    print(round(result, 3))

    ## bootstrap
    if (reps > 0) {
        Cvb = Cvw1b = Cvw2b = rep(NA, reps)
        for (b in 1:reps) {
            f.train = sample(1:n0, replace=TRUE)
            f.test  = sample(1:n1, replace=TRUE)

            Y1b   = Y1[f.test]
            mdlb  = update(mdl, data=dat.train[f.train,])
            scr1b = predict(mdlb, newdata=dat.test[f.test,])
            scr0b = mdlb$fitted.values

            xtrainb = xtrain[f.train,,drop=FALSE]
            xtestb  = xtest[f.test,,drop=FALSE]
            wgt1b   = densratio.appe(scr0b,   scr1b,  method, sigma, lambda, kernel_num, fold, stabilize, qstb)
            wgt2b   = densratio.appe(xtrainb, xtestb, method, sigma, lambda, kernel_num, fold, stabilize, qstb)

            Cvb[b]   = cvalest.bin(Y1b, scr1b)
            Cvw1b[b] = cvalest.bin(Y1b, scr1b, wgt1b)
            Cvw2b[b] = cvalest.bin(Y1b, scr1b, wgt2b)
        }

        ## se
        Cvse   = sd(Cvb,   na.rm=TRUE)
        Cvw1se = sd(Cvw1b, na.rm=TRUE)
        Cvw2se = sd(Cvw2b, na.rm=TRUE)

        ## percentile ci
        cl     = c((1-conf.level)/2, 1 - (1-conf.level)/2)
        Cvci   = quantile(Cvb,   cl, na.rm=TRUE)
        Cvw1ci = quantile(Cvw1b, cl, na.rm=TRUE)
        Cvw2ci = quantile(Cvw2b, cl, na.rm=TRUE)

        ## approx ci
        Cvcia   = Cv   + Cvse   * qnorm(cl)
        Cvw1cia = Cvw1 + Cvw1se * qnorm(cl)
        Cvw2cia = Cvw2 + Cvw2se * qnorm(cl)

        ## output
        message("\nPoint & Interval estimates:")
        result = cbind(result,
            c(Cvse, Cvw1se, Cvw2se),
            rbind(Cvci,  Cvw1ci,  Cvw2ci),
            rbind(Cvcia, Cvw1cia, Cvw2cia))
        names(result) = c('Est', 'SE', 'Percentile.l', 'Percentile.u', 'Approx.l', 'Approx.u')
        print(round(result, 3))
    }

    invisible(result)
}
