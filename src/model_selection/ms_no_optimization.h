//
// Created by jakob on 9/14/24.
//

#ifndef CPP_MS_NO_OPTIMIZATION_H
#define CPP_MS_NO_OPTIMIZATION_H

#include "ms_base.h"
#include "../mombf/modselIntegrals.h"

// Different model selection functions:
Eigen::VectorXi model_selection_no_optimization(
        Eigen::VectorXd y,
        Eigen::MatrixXi x,
        int n_iter,
        msPriorSpec prior_coef,
        msPriorSpec prior_delta,
        double phi,
        Eigen::VectorXi w_i,
        int n_observations,
        int n_timeperiods
);

#endif //CPP_MS_NO_OPTIMIZATION_H
