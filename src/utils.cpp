//
// Created by jakob on 4/13/25.
//
#include "utils.h"

namespace bisam {
    arma::ivec repeat_value(int value, int times) {
        arma::ivec result(times);
        result.fill(value);
        return result;
    }

    arma::mat kronecker_product(const arma::mat &A, const arma::mat &B) {
        arma::mat result(A.n_rows * B.n_rows, A.n_cols * B.n_cols);
        for (arma::uword i = 0; i < A.n_rows; ++i) {
            for (arma::uword j = 0; j < A.n_cols; ++j) {
                result.submat(i * B.n_rows, j * B.n_cols, (i + 1) * B.n_rows - 1, (j + 1) * B.n_cols - 1) = A(i, j) * B;
            }
        }
        return result;
    }
}
