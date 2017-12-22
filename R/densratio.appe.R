#' A wrapper function
#' A wrapper function to use "densratio" function from the densratio package.
#' @import stats
#' @importFrom densratio densratio
#' @param xtrain A dataframe used to construct a prediction model.
#' @param xtest A dataframe corresponding to a validation (testing) data.
#' @param method uLSIF or KLIEP. Same as the argument in densratio function from densratio package.
#' @param sigma A positive numeric vector corresponding to candidate values of a bandwidth for Gaussian kernel. Same as the argument in densratio function from densratio package.
#' @param lambda A positive numeric vector corresponding to candidate values of a regularization parameter. Same as the argument in densratio function from densratio package.
#' @param kernel_num A positive integer corresponding to number of kernels. Same as the argument in densratio function from densratio package.
#' @param fold A positive integer corresponding to a number of the folds of cross-validation in the KLIEP method. Same as the argument in densratio function from densratio package.
#' @param stabilize A logical value as to whether tail weight stabilization is performed or not. If TRUE, both tails of the estimated density ratio distribution are replaced by the constant value which is specified at qstb option.
#' @param qstb A positive numerical value less than 1 to control the degree of weight stabilization. Default value is 0.025, indicating estimated density ratio values less than the 2.5 percentile and more than the 97.5 percentile are set to 2.5 percentile and 97.5 percentile, respectively.
#' @return numerical vector
densratio.appe <-
function(xtrain, xtest, method="uLSIF", sigma=NULL, lambda=NULL, kernel_num=NULL, fold=5, stabilize=TRUE, qstb=0.025) {

    xtrain = as.matrix(xtrain)
    xtest  = as.matrix(xtest)

    if (is.null(kernel_num)) kernel_num = 100
    
    if (is.null(sigma)) {
        center = matrix(xtest[sample(1:nrow(xtest), kernel_num),], kernel_num, ncol(xtest))
        sigma  = as.array(quantile((dist(center))))
        sigma  = unique(sigma[ sigma>0.001 ])
    }

    if (is.null(lambda)) lambda = "auto"

    if (method == "uLSIF" || method == "KLIEP") {
        wgt = densratio(xtrain, xtest, method, sigma, lambda, kernel_num, fold, verbose=FALSE)$compute_density_ratio(xtest)
#    } else if (method == "gam") {
#        wgt = densratio.gam(xtrain, xtest, stabilize)
    } else {
#        stop("\n\nmethod should be either in ('uLSIF', 'KLIEP', 'gam').\n\n")
        stop("\n\nmethod should be either in ('uLSIF', 'KLIEP').\n\n")
    }

    ## tail-weight stabilization
    if (stabilize) {
        vl = quantile(wgt, qstb)
        wgt[ wgt < vl ] = vl
        vl = quantile(wgt, 1-qstb)
        wgt[ wgt > vl ] = vl
    }

    return(wgt)
}
