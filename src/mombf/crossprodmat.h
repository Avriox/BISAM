// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-

#ifndef CROSSPRODMAT
#define CROSSPRODMAT 1

#ifdef RCPP_EIGEN
// R environment
#include <RcppEigen.h>
// [[Rcpp::depends(RcppEigen)]]
#else
// Pure C++ environment
#include <Eigen/Dense>
#include <Eigen/Sparse>

#endif

// Standard library includes (common to both environments)
#include <memory>
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// #include "RcppArmadillo.h"

// via the depends attribute we tell Rcpp to create hooks for
// RcppArmadillo so that the build process will know what to do
//
// [[Rcpp::depends(RcppArmadillo)]]

// using namespace Rcpp;
using namespace std;

// CLASS crossprodmat using basic C++
// class crossprodmat {
//
// public:
//
//
//     crossprodmat(double *mymat, int nrowx, int ncolx, bool dense);
//
//     ~crossprodmat();
//
//     double at(int k);  //Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
//
//
// private:
//
//     double *x;           //X stored as a vector
//     int nrowx;
//     int ncolx;
//     int *userows; //optional slot indicating the indexes of the rows in x to be used when computing XtX. That is XtX= t(x[userows,]) %*% x[userows,]
//     int nuserows; //number of rows in x to be used when computing XtX
//     int userowsini; //if userows not provided, use x[userowsini : userowsini+nuserows-1,] to compute XtX (default userowsini=0)
//     bool dense; //if true then matrix is stored in XtXd, else in XtXs
//     double *XtXd;
//     arma::sp_mat XtXs;  //equivalent to SpMat<double> XtXs
//     arma::SpMat<short> XtXcomputed; //bool entries indicating if XtX has been computed
//
// };

class crossprodmat {
   public:
    crossprodmat(double *mymat, int nrowx, int ncolx, bool dense);

    ~crossprodmat();

    double at(int k);  // Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow

   private:
    //    std::unique_ptr<double[]> x;
    //    Eigen::MatrixXd XtXd;
    std::vector<double> XtXd;
    //    Eigen::SparseMatrix<double> XtXs;
    //    Eigen::SparseMatrix<short> XtXcomputed;
};

#endif
