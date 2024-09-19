#include "crossprodmat.h"

using namespace std;
//using namespace arma;
//
//
//crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense) {
//
//    this->nrowx = nrowx;
//    this->ncolx = ncolx;
//    this->userowsini = 0;
//    this->nuserows = nrowx;
//    this->userows = NULL;
//
//    if (dense) {
//        this->XtXd = mymat;
//        this->dense = true;
//    } else {
//        this->x = mymat;
//        this->dense = false;
//        (this->XtXs) = arma::sp_mat(ncolx, ncolx);
//        (this->XtXcomputed) = arma::SpMat<short>(ncolx, ncolx);
//    }
//}
//
//
////Class destructor
//crossprodmat::~crossprodmat() {}
//
//
////Access element with vector-type index A(k)= A(i,j) where j= k/nrow; i= k % nrow
//double crossprodmat::at(int k) {
//
//    if (dense) {
//
//        return XtXd[k];
//
//    } else {
//
//        int i = k % ncolx, j = k / ncolx;
//
//        if (XtXcomputed.at(i, j) == 0) {  //if this entry has not been already computed
//
//            int iini, jini, k;
//            double ans = 0;
//            if (this->userows == NULL) {
//
//                for (k = this->userowsini, iini = i * nrowx, jini = j * nrowx;
//                     k < this->nuserows + this->userowsini; k++)
//                    ans += x[k + iini] * x[k + jini];
//
//            } else {
//
//                for (k = 0, iini = i * nrowx, jini = j * nrowx; k < this->nuserows; k++)
//                    ans += x[(this->userows[k]) + iini] * x[(this->userows[k]) + jini];
//
//            }
//
//            XtXcomputed(i, j) = 1;
//            XtXs(i, j) = ans;
//        }
//
//        return XtXs.at(i, j);
//
//    }
//
//
//}

#include "crossprodmat.h"
#include <stdexcept>
#include <algorithm>

crossprodmat::crossprodmat(double *mymat, int nrowx, int ncolx, bool dense)
        : XtXd(mymat, mymat + ncolx * ncolx) {

//    if (dense) {
//        XtXd = Eigen::Map<Eigen::MatrixXd>(mymat, ncolx, ncolx);
//    } else {
//        x = std::make_unique<double[]>(nrowx * ncolx);
//        std::copy(mymat, mymat + (nrowx * ncolx), x.get());
//        XtXs.resize(ncolx, ncolx);
//        XtXcomputed.resize(ncolx, ncolx);
//    }
}

crossprodmat::~crossprodmat() = default;

// Original Version
//double crossprodmat::at(int k) {
//    if (k < 0 || k >= ncolx * ncolx) {
//        throw std::out_of_range("Index out of bounds");
//    }
//
//    if (dense) {
//        int i = k % ncolx, j = k / ncolx;
//        return XtXd(i, j);
//    } else {
//        int i = k % ncolx, j = k / ncolx;
//
//        if (XtXcomputed.coeff(i, j) == 0) {  // if this entry has not been already computed
//            double ans = 0.0;
//            if (userows == nullptr) {
//                for (int k = userowsini; k < nuserows + userowsini; ++k) {
//                    ans += x[k + i * nrowx] * x[k + j * nrowx];
//                }
//            } else {
//                for (int k = 0; k < nuserows; ++k) {
//                    ans += x[userows[k] + i * nrowx] * x[userows[k] + j * nrowx];
//                }
//            }
//
//            XtXcomputed.coeffRef(i, j) = 1;
//            XtXs.coeffRef(i, j) = ans;
//        }
//
//        return XtXs.coeff(i, j);
//    }
//}

double crossprodmat::at(int k) {
//    if (k < 0 || k >= ncolx * ncolx) {
//        throw std::out_of_range("Index out of bounds");
//    }

//    int i = k % ncolx;
//    int j = k / ncolx;
//    return XtXd(i, j);

    return XtXd[k];

}

