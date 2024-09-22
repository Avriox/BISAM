//
// Created by jakob on 6/25/24.
//

#include "b_ism.h"
#include "exec_timer.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

void b_ism(
        const Eigen::MatrixXd &data,
        bool include_constant,
        bool tfe,
        bool ife,
        bool iis,
        bool sis,
        int y_index,
        int i_index,
        int t_index,
        long Ndraw,
        long Nburn,
        std::string b_prior,
        double lambda_b,
        double c0,
        double C0,
        bool geweke,
        BismResults &results,
        ModelSelectionVersion model_selection_version
) {

//    if (geweke) {
//        lambda_b = lambda_g = 1.234;
//        Ndraw = 100000;
//        c0 = 3;
//        C0 = 1;
//        sis = true;
//        b_prior = "g";
//    }


    // Create column names
    // TODO need to check whether we even need this
//    std::vector<std::string> colnames = {"i", "t", "y"};
//    for (int i = 0; i < colnames.size(); ++i) {
//        if (i != y_index && i != i_index && i != t_index) {
//            colnames[i] = "beta." + colnames[i];
//        }
//    }

    FunctionTimer timer;
    timer.start_section("Preparation");

    // Create a vector of column indices to keep
    // Keep all of the columns that are not y, i or t_index
    std::vector<int> indices_to_keep;
    for (int i = 0; i < data.cols(); ++i) {
        if (i != y_index && i != i_index && i != t_index) {
            indices_to_keep.push_back(i);
        }
    }

    // Claud said y_index - 1 but that is wrong!
    Eigen::MatrixXd y = data.col(y_index);

    // Extract X and X_
    Eigen::MatrixXd X(data.rows(), indices_to_keep.size());
    for (int i = 0; i < indices_to_keep.size(); ++i) {
        X.col(i) = data.col(indices_to_keep[i]);
    }


    // This is a copy as it should be!
//    Eigen::MatrixXd X_ = X;

    // Find unique values for n and t
    std::set<double> unique_i;
    std::set<double> unique_t;

    for (int i = 0; i < data.rows(); ++i) {
        unique_i.insert(data(i, i_index));
        unique_t.insert(data(i, t_index));
    }

    int n = unique_i.size();
    int t = unique_t.size();
    int N = n * t;

    // Build X matrix
    if (include_constant) {
        Eigen::MatrixXd const_col = Eigen::MatrixXd::Ones(X.rows(), 1);
        X.conservativeResize(X.rows(), X.cols() + 1);
        X.col(X.cols() - 1) = const_col;
    }

//     TODO UNTESTED
    if (ife) {
        Eigen::MatrixXd IFE = kroneckerProduct(Eigen::MatrixXd::Identity(n, n), Eigen::MatrixXd::Ones(t, 1));
        X.conservativeResize(X.rows(), X.cols() + IFE.cols());
        X.rightCols(IFE.cols()) = IFE;
    }
    // TODO UNTESTED
    if (tfe) {
        Eigen::MatrixXd TFE = kroneckerProduct(Eigen::MatrixXd::Ones(n, 1), Eigen::MatrixXd::Identity(t, t));
        X.conservativeResize(X.rows(), X.cols() + TFE.cols());
        X.rightCols(TFE.cols()) = TFE;
    }
    // TODO UNTESTED
    if (tfe && ife && !include_constant) {
        std::cerr
                << "Warning: Both time and unit fixed effects used.\nDropping first indiv. FE to avoid perfect collinearity"
                << std::endl;
        X.conservativeResize(X.rows(), X.cols() - 1);
    }

    // Find Matrix dimensions
    int p = X.cols();
    int r = n * (t - 1 - iis);

    // Creating z and Z needs to be done in multiple steps here
    // z <- lower.tri(matrix(1,nrow=t,ncol=t),diag = T)[,-c(1,t*iis)]
    // We first create a "perfect" lower triangular matrix and then we make sure the
    // first row has only zeros and the last column ends with 2 ones
    Eigen::MatrixXi z = Eigen::MatrixXi::Zero(t, t);

    for (int i = 0; i < t; ++i) {
        for (int j = 0; j <= i; ++j) {
            z(i, j) = 1;
        }
    }



    // Remove the first column and the column corresponding to t*iis
    // TODO validate that -1 here is correct (zero based indexing)
    // TODO According to claude this may need to be changed to t * iis. Not tested yet
    // TODO I think it should be correct as is, looking at the ouput.
    int col_to_remove1 = 0;
    int col_to_remove2 = t * iis - 1;

    Eigen::MatrixXi z_reduced(t, t - 2);

    int col_idx = 0;
    for (int j = 0; j < t; ++j) {
        if (j != col_to_remove1 && j != col_to_remove2) {
            z_reduced.col(col_idx++) = z.col(j);
        }
    }

    z = z_reduced;

    // Kronecker product with identity matrix
    Eigen::MatrixXi I = Eigen::MatrixXi::Identity(n, n);
    Eigen::MatrixXi Z = Eigen::MatrixXi::Zero(n * z.rows(), n * z.cols());

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (I(i, j) == 1) {
                Z.block(i * z.rows(), j * z.cols(), z.rows(), z.cols()) = z;
            }
        }
    }



    // This also needs to be broken up significantly
    // sis_grid <- expand.grid(unique(data[, t_index])[-c(1,t*iis)], unique(data[, i_index]))
    // unique values have already been calculated

    // Generate unique values for the specified columns
    std::vector<double> unique_t_vector(unique_t.begin(), unique_t.end());
    std::vector<double> unique_i_vector(unique_i.begin(), unique_i.end());

    // Remove the first element and the element at t*iis - 1 from unique_t_vector
    if (t * iis - 1 < unique_t_vector.size()) {
        unique_t_vector.erase(unique_t_vector.begin() + (t * iis - 1));
    }
    unique_t_vector.erase(unique_t_vector.begin());

    // Create the grid (combination of unique_t and unique_i)
    std::vector<std::pair<double, double>> sis_grid;
    for (double i: unique_i_vector) {
        for (double t: unique_t_vector) {
            sis_grid.emplace_back(t, i);
        }
    }


    // Create column names for Z
    // TODO probably not useful right now but maybe when we send data back to R
//    std::vector<std::string> colnames_Z;
//    for (const auto &pair: sis_grid) {
//        colnames_Z.push_back(
//                "sis." + std::to_string((int) pair.first) + "." + std::to_string((int) pair.second));
//    }


    // Perform the cross-product
    // TODO do we need to check whether XX is positive for the decomp to work?
    // TODO I think it is always pos def
    Eigen::MatrixXd XX = X.transpose() * X;
    Eigen::MatrixXi ZZ = Z.transpose() * Z;

//    std::cout << "X: \n" << X << std::endl << std::endl;
//
//    std::cout << "XX: \n" << XX << std::endl << std::endl;

    // Prior Parameters
    Eigen::VectorXd b0;
//    Eigen::MatrixXd B0(p, p);
//    Eigen::MatrixXd B0_inv(p, p);
    Eigen::LDLT<Eigen::MatrixXd> ldltOfXX(XX);


//    if (b_prior == "g") {
    b0 = Eigen::VectorXd::Zero(p);
//    } else if (b_prior == "f") {
//     TODO UNTESTED BRANCH NEVER REACHED IN TESTING
//        b0 = ldltOfXX.solve(X.transpose() * y);
//        B0 = ldltOfXX.solve(Eigen::MatrixXd::Identity(p, p)) * lambda_b;
//    }

//    B0 = ldltOfXX.solve(Eigen::MatrixXd::Identity(p, p)) * lambda_b;
//    B0_inv = XX / lambda_b;

    Eigen::VectorXd g0 = Eigen::VectorXd::Zero(r);

    // Store
    long Nstore;
    if (Nburn >= Ndraw) {
        throw std::runtime_error("The number of burn-in exceeds number of draws");
    } else {
        Nstore = Ndraw - Nburn;
    }

    // Main
    Eigen::MatrixXd b_store = Eigen::MatrixXd::Constant(Nstore, p, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd g_store = Eigen::MatrixXd::Constant(Nstore, r, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd s2_store = Eigen::MatrixXd::Constant(Nstore, 1, std::numeric_limits<double>::quiet_NaN());

    // Setting col names is not relly a thing here!
    // TODO check whether I need to set the col names for:
    //  b_store, g_store, s2_store, w_store, w_pip_store
    //  before returning them to R

    // selection
    Eigen::MatrixXi w_store = Eigen::MatrixXi::Constant(Nstore, r, std::numeric_limits<int>::quiet_NaN());

    // hyper
    Eigen::MatrixXd v_store = Eigen::MatrixXd::Constant(Nstore, 1, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd w_pip_store = Eigen::MatrixXd::Constant(Nstore, r, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd xi_store = Eigen::MatrixXd::Constant(Nstore, 1, std::numeric_limits<double>::quiet_NaN());
    Eigen::MatrixXd tau_store = Eigen::MatrixXd::Constant(Nstore, 1, std::numeric_limits<double>::quiet_NaN());

    // G0_store is a 3-dimensional array
    // TODO this matrix might need checking!
    std::vector<Eigen::MatrixXd> G0_store(Nstore,
                                          Eigen::MatrixXd::Constant(r, r, std::numeric_limits<double>::quiet_NaN()));

    // Starting Values
    Eigen::VectorXd b_i = Eigen::VectorXd::Ones(p);
    Eigen::VectorXd g_i = Eigen::VectorXd::Ones(r); // r
    Eigen::VectorXd g_i_n0 = g_i; // TODO besserer name
    double s2_i = 1.0;

    // TODO make int
    Eigen::VectorXi w_i = Eigen::VectorXi::Zero(r);

    // creating z_cols is a bit more complicated
    Eigen::VectorXi z_cols(n * (t - 1 - iis));
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < (t - 1 - iis); ++j) {
            z_cols(i * (t - 1 - iis) + j) = j + 1;
        }
    }

//    Eigen::VectorXd w_1 = z_cols.cast<double>().array() * w_i.array(); // columns where w_i is 1


    double tau_g = 1.0;
    double xi_g = 1.0;
    Eigen::VectorXd lamb_g = Eigen::VectorXd::Ones(r);


    // Prep Calculations
    double cN = c0 + N / 2.0;
    // TODO didnt we calculate that already?
    Eigen::MatrixXd XX_inv = XX.inverse();
    Eigen::VectorXd Pr_g = Eigen::VectorXd::Constant(n * (t - 1 - iis), std::numeric_limits<double>::quiet_NaN());


    // Proximity Matrix for Conditional break probability
    // TODO move / make parameter?
    int cutoff = 4; // number of lags allowed to matter
    double prob = 0.8; // probability of after break being no break

    Eigen::MatrixXd NNMat = Eigen::MatrixXd::Zero(r / n, r / n);
    for (int i = 0; i < r / n; ++i) {
        for (int j = 0; j < r / n; ++j) {
            if (i != j) {
                int distance = std::abs(i - j);
                if (distance <= cutoff) {
                    NNMat(i, j) = std::pow(1 - prob, distance - 1);
                }
            }
        }
    }


    // Set which version of Model Selection to use
    // All of them HAVE TO take the same parameters and return the same type!
    // Using std::function to store the selected function
    std::function<Eigen::VectorXi(Eigen::VectorXd, Eigen::MatrixXi, int, msPriorSpec, msPriorSpec,
                                  double, Eigen::VectorXi, int n_observations, int n_timeperiods)> model_selection;


    switch (model_selection_version) {
        case ModelSelectionVersion::NO_OPTIMIZATION: {
            model_selection = model_selection_no_optimization;
            break;
        };
        case ModelSelectionVersion::SPLIT_MATRIX: {
            model_selection = model_selection_split_z;
            break;
        };
        case ModelSelectionVersion::SPLIT_MATRIX_PARALLEL: {
            model_selection = model_selection_parallel_z;
        }
    }

    // Create a progress bar for R
    // Create progress bar
    progressbar bar(Ndraw);
    bar.set_todo_char(" ");
    bar.set_done_char("■");
    bar.set_opening_bracket_char("[");
    bar.set_closing_bracket_char("]");


    // Initialize random number generator
    std::random_device rd;
    std::mt19937 rng(rd());

    // Allocate the answer vector for rnlp
    // TODO can probably be allocated once outside of the loop
    Eigen::VectorXd ans = Eigen::VectorXd::Zero(Z.cols());

    timer.end_section("Preparation");

    // With default values i will go from -999 to 9000
    for (int i = 1 - Nburn; i <= Nstore; ++i) {

        //================ Geweke Test =================.====
        // Für Jakob: witziger test um richtigkeit des Gibbs-Samplers zu überprüfen
        // TODO UNTESTED!
//        if (geweke) {
//            Eigen::VectorXd e = Eigen::VectorXd::NullaryExpr(N, [&]() {
//                return std::normal_distribution<>(0, std::sqrt(s2_i))(rng);
//            });
//            y = X * b_i + Z.cast<double>() * g_i + e;
//        }

        //=================== draw p(s2|a,b,g,y) ==========.====
        // samplen der varianz gegeben der andere
        // n
        // cN

        timer.start_section("Sample S2");

        // TODO remove casting?
        double CN = C0 + 0.5 * (y - X * b_i - Z.cast<double>() * g_i).squaredNorm();
        s2_i = 1.0 / std::gamma_distribution<>(cN, 1.0 / CN)(rng);

        timer.end_section("Sample S2");
        //=================== draw p(b|a,g,s2,y) ==========.====
        // samplen der betas

        timer.start_section("Sample b");
        Eigen::MatrixXd BN_inv = (1 / s2_i + 1 / lambda_b) * XX;
        Eigen::MatrixXd BN = s2_i * lambda_b / (s2_i + lambda_b) * XX_inv;
        Eigen::VectorXd bN = 1 / (s2_i + lambda_b) *
                             (lambda_b * XX_inv * (X.transpose() * (y - Z.cast<double>() * g_i)) + s2_i * b0);


        // TODO I am assuming his block to be correct, the results look very comparable to the results I get in the R code
        Eigen::VectorXd b_i_try;
//        try {
        Eigen::LLT<Eigen::MatrixXd> cholBN(BN);
        if (cholBN.info() != Eigen::Success) {
            throw std::runtime_error("Cholesky decomposition of BN failed");
        }
        Eigen::VectorXd z = Eigen::VectorXd::NullaryExpr(p,
                                                         [&]() {
                                                             return std::normal_distribution<>(0, 1)(rng);
                                                         });
        b_i_try = bN + cholBN.matrixU().transpose() * z;
//        } catch (const std::exception &e) {
//            // TODO the fallback has not been validated
//            // Fallback method: use multivariate normal distribution directly
//            std::normal_distribution<> normal(0, 1);
//            Eigen::VectorXd z(p);
//            for (int i = 0; i < p; ++i) {
//                z(i) = normal(rng);
//            }
//            b_i_try = bN + BN.llt().matrixL() * z;
//        }
        b_i = b_i_try;

        timer.end_section("Sample b");

        // Prior Variances:
        // TODO brauchen wir G0 überhaupt?
        Eigen::VectorXd diag_elements = (lamb_g.array() * tau_g * s2_i).matrix();
        Eigen::MatrixXd G0 = diag_elements.asDiagonal();


//        if (sis) {

        timer.start_section("Model selection");

        Eigen::VectorXd y_hat = y - X * b_i;

        // ausm loop
        msPriorSpec priorCoef = imomprior(tau_g);
        msPriorSpec priorDelta = modelbbprior(1, 1);


//        std::cout << "Z:" << std::endl << Z << std::endl << std::endl;
//
//        std::cout << "w_i:" << std::endl << w_i << std::endl << std::endl;


        // Use the correct version of model selection as defined earlier
        // TODO I am not 100% sure that col_idx is correct for n_observations
        w_i = model_selection(y_hat, Z, 2, priorCoef, priorDelta, s2_i, w_i, col_idx, t);

        timer.end_section("Model selection");


        timer.start_section("RNLP");

        // Prepare a Z that only has the cols where w_i is non zero
        // Count non-zero elements in w_i
        int nonZeroCount = (w_i.array() != 0).count();

        // Create a new matrix with the same number of rows as Z and columns equal to nonZeroCount
        Eigen::MatrixXi filteredZ(Z.rows(), nonZeroCount);

        Eigen::VectorXd filteredGi(nonZeroCount);

        int currentCol = 0;
        for (int i = 0; i < w_i.size(); ++i) {
            if (w_i(i) != 0) {
                // Copy the corresponding column from Z to filteredZ
                filteredZ.col(currentCol) = Z.col(i);
                filteredGi[currentCol] = g_i_n0(i);
                ++currentCol;
            }
        }

        rnlpPost_lm(ans.data(), 1, 0, 1, y_hat.data(), filteredZ.data(), y_hat.size(), filteredZ.cols(),
                    1, tau_g, c0, C0, 1, filteredGi.data(), s2_i);

        g_i = Eigen::VectorXd::Zero(r);

        int index = 0;
        for (int h = 0; h < w_i.size(); ++h) {
            if (w_i(h) != 0) {
                g_i(h) = ans(index);
                g_i_n0(h) = ans(index);
                index++;
            }
        }

        timer.end_section("RNLP");

//        }


        // Store results
        if (i > 0) {
            timer.start_section("Store");
            // Regular
            b_store.row(i - 1) = b_i.transpose();
            g_store.row(i - 1) = g_i.transpose();
            w_store.row(i - 1) = w_i.transpose();
            s2_store(i - 1) = s2_i;

            // Hyper
//                v_store(i - 1) = v_i;
            w_pip_store.row(i - 1) = Pr_g.transpose();
            xi_store(i - 1) = xi_g;
            tau_store(i - 1) = tau_g;

            // Store G0 matrix
            for (int j = 0; j < r; ++j) {
                for (int k = 0; k < r; ++k) {
                    G0_store[i - 1](j, k) = G0(j, k);
                }
            }

            timer.end_section("Store");
        }


        // Update the R progress bar
        bar.update();
    }



    // Populate the results
    results.b_store = b_store;
    results.s2_store = s2_store;
    results.g_store = g_store;
    results.w_store = w_store;
    results.v_store = v_store;
    results.w_pip_store = w_pip_store;


    timer.print_section_summary();
}


// Wrapper function for R
// [[Rcpp::export]]
Rcpp::List b_ism_wrapper(
        Rcpp::NumericMatrix data,
        bool include_constant,
        bool tfe,
        bool ife,
        bool iis,
        bool sis,
        int y_index,
        int i_index,
        int t_index,
        long Ndraw,
        long Nburn,
        std::string b_prior,
        double lambda_b,
        double c0,
        double C0,
        bool geweke,
        std::string model_sel_optimization
) {
    try {
        // Convert R's NumericMatrix to Eigen::MatrixXd
        Eigen::MatrixXd data_eigen = Rcpp::as<Eigen::MatrixXd>(data);

        // Create a BismResults structure to store results
        BismResults results;

        ModelSelectionVersion model_selection_version;
        if (model_sel_optimization == "No Optimization") {
            std::cout << "Running with model_selection: 'No Optimization'" << std::endl;
            model_selection_version = ModelSelectionVersion::NO_OPTIMIZATION;

        } else if (model_sel_optimization == "Split Z") {
            model_selection_version = ModelSelectionVersion::SPLIT_MATRIX;
            std::cout << "Running with model_selection: 'Split Z Matrix'" << std::endl;

        } else if (model_sel_optimization == "Parallel Z") {
            model_selection_version = ModelSelectionVersion::SPLIT_MATRIX_PARALLEL;
            std::cout << "Running with model_selection: 'Split and Parallel Z Matrix'" << std::endl;

        } else {
            std::cout << "Selected Model Selection Version does not exist. Defaulting to 'No Optimization'";
            model_selection_version = ModelSelectionVersion::NO_OPTIMIZATION;
        }



        // Call the b_ism function
        b_ism(
                data_eigen, include_constant, tfe, ife, iis, sis, y_index, i_index, t_index,
                Ndraw, Nburn, b_prior, lambda_b, c0, C0, geweke, results, model_selection_version
        );

        // Convert BismResults to Rcpp::List
        Rcpp::List out = Rcpp::List::create(
                Rcpp::Named("draws") = Rcpp::List::create(
                        Rcpp::Named("beta") = results.b_store,
                        Rcpp::Named("sigma2") = results.s2_store,
                        Rcpp::Named("sis") = results.g_store,
                        Rcpp::Named("omega") = results.w_store
                ),
                Rcpp::Named("meta") = Rcpp::List::create(
                        Rcpp::Named("MCMC") = Rcpp::NumericVector::create(
                                Rcpp::Named("Ndraw") = Ndraw,
                                Rcpp::Named("Nburn") = Nburn,
                                Rcpp::Named("Nstore") = Ndraw - Nburn
                        )
                )
        );

        if (iis && sis) {
            out["hyper_draws"] = Rcpp::List::create(
                    Rcpp::Named("nu") = results.v_store,
                    Rcpp::Named("omega_pip") = results.w_pip_store
            );
        }

        return out;
    } catch (const std::exception &ex) {
        Rcpp::stop("Exception caught: %s", ex.what());
    } catch (...) {
        Rcpp::stop("Unknown exception caught");
    }
}