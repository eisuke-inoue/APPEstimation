#' L1 and L2 errors adjusted for predictor distributions
#' Calculates adjusted L1 and L2 errors by predictor distributions for a linear model.
#' @export
#' @import stats
#' @param mdl A lm object describing a prediction model to be evaluated.
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
appe.lm <-
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
    L1   = mean(abs(Y1 - scr1))
    L1w1 = weighted.mean(abs(Y1 - scr1), w=wgt1)
    L1w2 = weighted.mean(abs(Y1 - scr1), w=wgt2)
    L2   = mean((Y1 - scr1)^2)
    L2w1 = weighted.mean((Y1 - scr1)^2, w=wgt1)
    L2w2 = weighted.mean((Y1 - scr1)^2, w=wgt2)

    message("\nPoint estimates:")
    result = data.frame(c(L1, L1w1, L1w2, L2, L2w1, L2w2))
    names(result) = 'Estimate'
    row.names(result) = c('L1','L1 adjusted by score','L1 adjusted by predictors','L2','L2 adjusted by score','L2 adjusted by predictors')
    print(round(result, 3))

    ## bootstrap
    if (reps > 0) {
        L1b = L1w1b = L1w2b = L2b = L2w1b = L2w2b = rep(NA, reps)
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

            L1b[b]   = mean(abs(Y1b - scr1b))
            L1w1b[b] = weighted.mean(abs(Y1b - scr1b), w=wgt1b)
            L1w2b[b] = weighted.mean(abs(Y1b - scr1b), w=wgt2b)

            L2b[b]   = mean((Y1b - scr1b)^2)
            L2w1b[b] = weighted.mean((Y1b - scr1b)^2, w=wgt1b)
            L2w2b[b] = weighted.mean((Y1b - scr1b)^2, w=wgt2b)
        }

        ## se
        L1se   = sd(L1b,   na.rm=TRUE)
        L1w1se = sd(L1w1b, na.rm=TRUE)
        L1w2se = sd(L1w2b, na.rm=TRUE)
        L2se   = sd(L2b,   na.rm=TRUE)
        L2w1se = sd(L2w1b, na.rm=TRUE)
        L2w2se = sd(L2w2b, na.rm=TRUE)

        ## percentile ci
        cl     = c((1-conf.level)/2, 1 - (1-conf.level)/2)
        L1ci   = quantile(L1b,   cl, na.rm=TRUE)
        L1w1ci = quantile(L1w1b, cl, na.rm=TRUE)
        L1w2ci = quantile(L1w2b, cl, na.rm=TRUE)
        L2ci   = quantile(L2b,   cl, na.rm=TRUE)
        L2w1ci = quantile(L2w1b, cl, na.rm=TRUE)
        L2w2ci = quantile(L2w2b, cl, na.rm=TRUE)

        ## approx ci
        L1cia   = L1   + L1se   * qnorm(cl)
        L1w1cia = L1w1 + L1w1se * qnorm(cl)
        L1w2cia = L1w2 + L1w2se * qnorm(cl)
        L2cia   = L2   + L2se   * qnorm(cl)
        L2w1cia = L2w1 + L2w1se * qnorm(cl)
        L2w2cia = L2w2 + L2w2se * qnorm(cl)

        ## output
        message("\nPoint & Interval estimates:")
        result = cbind(result,
            c(L1se, L1w1se, L1w2se, L2se, L2w1se, L2w2se),
            rbind(L1ci,  L1w1ci,  L1w2ci,  L2ci,  L2w1ci,  L2w2ci),
            rbind(L1cia, L1w1cia, L1w2cia, L2cia, L2w1cia, L2w2cia))
        names(result) = c('Estimate', 'Std.Error', 'Percentile.l', 'Percentile.u', 'Approx.l', 'Approx.u')
        print(round(result, 3))
    }

    invisible(result)
}
