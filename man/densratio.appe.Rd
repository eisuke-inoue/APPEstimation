% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/densratio.appe.R
\name{densratio.appe}
\alias{densratio.appe}
\title{A wrapper function
A wrapper function to use "densratio" function from the densratio package.}
\usage{
densratio.appe(xtrain, xtest, method = "uLSIF", sigma = NULL,
  lambda = NULL, kernel_num = NULL, fold = 5, stabilize = TRUE,
  qstb = 0.025)
}
\arguments{
\item{xtrain}{A dataframe used to construct a prediction model.}

\item{xtest}{A dataframe corresponding to a validation (testing) data.}

\item{method}{uLSIF or KLIEP. Same as the argument in densratio function from densratio package.}

\item{sigma}{A positive numeric vector corresponding to candidate values of a bandwidth for Gaussian kernel. Same as the argument in densratio function from densratio package.}

\item{lambda}{A positive numeric vector corresponding to candidate values of a regularization parameter. Same as the argument in densratio function from densratio package.}

\item{kernel_num}{A positive integer corresponding to number of kernels. Same as the argument in densratio function from densratio package.}

\item{fold}{A positive integer corresponding to a number of the folds of cross-validation in the KLIEP method. Same as the argument in densratio function from densratio package.}

\item{stabilize}{A logical value as to whether tail weight stabilization is performed or not. If TRUE, both tails of the estimated density ratio distribution are replaced by the constant value which is specified at qstb option.}

\item{qstb}{A positive numerical value less than 1 to control the degree of weight stabilization. Default value is 0.025, indicating estimated density ratio values less than the 2.5 percentile and more than the 97.5 percentile are set to 2.5 percentile and 97.5 percentile, respectively.}
}
\value{
numerical vector
}
\description{
A wrapper function
A wrapper function to use "densratio" function from the densratio package.
}
