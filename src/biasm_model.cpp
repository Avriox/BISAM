//
// Created by jakob on 4/13/25.
//
#include "biasm_model.h"

#include "modelselection_strategy.h"

namespace bisam {
    // TODO bisam_regression
    //  * include_indiv_fixed_eff
    //  * improve variable names
    //  * use_phiinit -> sample phi extra
    //  * lambda_b -> prio_coef_var
    //
    // Store values in model selection:
    //


    BisamResult estimate_model(arma::mat &data,
                               int i_index,
                               int t_index,
                               int y_index,
                               long long num_draws,
                               long long num_burnin,
                               std::string b_prior,
                               double lambda_b,
                               double c0,
                               double C0,
                               double va,
                               double vb,
                               double tau,
                               bool use_phiinit,
                               bool const_val,
                               bool ife,
                               bool tfe,
                               bool iis,
                               bool sis,
                               int new_par_method,
                               int new_par_hesstype,
                               int new_par_optim_method,
                               int new_par_optim_maxit,
                               int new_par_B,
                               int new_par_knownphi,
                               int new_par_r,
                               double new_par_alpha,
                               double new_par_lambda,
                               ComputationStrategy strategy) {
        // [DEV]
        // This sets up a "random device" for development purposes. We need it because RcppArmadillo tries to
        // use R's random number generator per default which we have not loaded when we are developing in cpp.
#ifdef DEV

        std::random_device rd;
        std::mt19937 cpp_rng_generator(rd());
        std::normal_distribution<> cpp_normal_dist(0.0, 1.0); // Mean 0, StdDev 1
#endif


        /* -------------------------------------------------------------------------- */
        /*                               Variable Setup                               */
        /* -------------------------------------------------------------------------- */

        // This section does general setup of variables / settings - this is only needed
        // for development to test the code. Usually these values should be provided from
        // within R

        /* --------------- Setup for the remaining / gibbs parameters --------------- */
        // TODO descriptions are missing

        // int i_index         = 0; // Adjusting to 0-based indexing in C++
        // int t_index         = 1; // Adjusting to 0-based indexing in C++
        // int y_index         = 2; // Adjusting to 0-based indexing in C++
        // long long Ndraw     = 5000;
        // long long Nburn     = 500;
        // std::string b_prior = "g";
        // double lambda_b     = 100.0;
        // double c0           = 0.001;
        // double C0           = 0.001;
        // double va           = 1.0;
        // double vb           = 1.0;
        // double tau          = 1.0;
        // bool geweke         = false;
        // bool use_phiinit    = true;


        /* -------------------------------------------------------------------------- */
        /*                              Data Preparation                              */
        /* -------------------------------------------------------------------------- */

        // Prepare the generated data to get it into a form we need

        /* ---------------------------- Data Extraction ----------------------------- */

        // Extract the y vector from the sim data
        arma::vec y = data.col(y_index);

        // Extract initial X (excluding y, individual index, and time index)
        // Vector of which columns to exclude {0,1,2}
        arma::ivec exclude_cols = {
            y_index, i_index, t_index
        };

        // Our X matrix should have cols(data) - length(excluded_cols) columns
        int num_include_cols = data.n_cols - exclude_cols.n_elem;

        // Pre allocate the Matrix
        arma::mat X_(data.n_rows, num_include_cols);

        // Helper to keep track of the current insertion index into X
        int current_x_col = 0;

        // Loop over all columns of data. If the columns index is one of the ones in
        // exclude_cols, skip, else take the column and insert it into X at current col
        // index. Increment index.
        for (int j = 0; j < data.n_cols; ++j) {
            bool exclude = false;
            for (int k = 0; k < exclude_cols.n_elem; ++k) {
                if (j == exclude_cols(k)) {
                    exclude = true;
                    break;
                }
            }
            if (!exclude) {
                X_.col(current_x_col) = data.col(j);
                current_x_col++;
            }
        }

        // TODO we probably do not need a copy here???
        arma::mat X = X_; // Create a copy for further modifications

        // Determine n and t (number of units and time periods)
        arma::vec unique_n_vals = arma::unique(data.col(i_index));
        arma::vec unique_t_vals = arma::unique(data.col(t_index));
        int n                   = unique_n_vals.n_elem;
        int t                   = unique_t_vals.n_elem;
        long N                  = (n) * t;


        /* ----------------------------- Build X Matrix ----------------------------- */
        // Add columns depending on what to include.

        if (const_val) {
            X.insert_cols(X.n_cols, arma::ones(data.n_rows));
        }

        if (ife) {
            // Individual fixed effects
            arma::mat IFE = kronecker_product(arma::eye(n, n), arma::ones(t, 1));
            // No direct equivalent to R's colnames assignment here, but we can keep track of the meaning
            if (const_val) {
                IFE = IFE.cols(1, IFE.n_cols - 1);
            }
            X.insert_cols(X.n_cols, IFE);
        }

        if (tfe) {
            // Time fixed effects
            arma::mat TFE = kronecker_product(arma::ones(n, 1), arma::eye(t, t));
            // No direct equivalent to R's colnames assignment here
            if (const_val) {
                TFE = TFE.cols(1, TFE.n_cols - 1);
            }
            X.insert_cols(X.n_cols, TFE);
        }

        if (tfe && ife && !const_val) {
            std::cerr <<
                    "Warning: Both time and unit fixed effects used.\nDropping first indiv. FE to avoid perfect colinearity"
                    << std::endl;
            if (!X_.is_empty()) {
                X.shed_col(X_.n_cols); // Assuming the IFE dummies were added after the original X_
            } else if (const_val) {
                // If no original regressors, and IFE was added, it would be at the beginning
                X.shed_col(0);
            } else if (ife) {
                // If only IFE was added, it would be at the beginning
                X.shed_col(0);
            }
        }


        /* ------------------------- Find Matrix dimensions ------------------------- */

        int p = X.n_cols; // Find matrix dimensions

        // Construct z matrix
        arma::mat z_temp(t, t, arma::fill::ones);
        arma::mat z_lower_tri = arma::trimatl(z_temp);
        arma::mat z           = z_lower_tri.cols(2, t - 2);
        // Remove first two and the last column (adjusting indices for 0-based)
        // z = arma::conv_to<arma::mat>::from(arma::round(z));
        // Convert to integer (Armadillo doesn't have direct integer matrix type for mat)


        // Construct Z matrix
        // TODO I think Arma has no integer matrix type? Does that have performance implications?
        arma::mat Z = kronecker_product(arma::eye(n, n), z);
        Z           = arma::conv_to<arma::mat>::from(arma::round(Z)); // Convert to integer


        int r        = Z.n_cols;
        arma::mat XX = arma::trans(X) * X;
        arma::mat ZZ = arma::trans(Z) * Z;

        /* ---------------------------- Prior Parameters ---------------------------- */

        arma::vec b0;
        if (b_prior == "g") {
            b0.zeros(p);
        } else if (b_prior == "f") {
            // b0 = arma::inv(XX) * X.t() * y;
            // Note: For numerical stability, Armadillo provides a more direct way to solve this system:
            b0 = arma::solve(XX, arma::trans(X) * y);
        }

        arma::mat B0     = XX.i() * lambda_b;
        arma::mat B0_inv = XX / lambda_b;


        /* ---------------------------- Starting Values ----------------------------- */

        // arma::vec b_i(p, arma::fill::ones);
        // TODO reuse value from before
        arma::vec b_i = arma::solve(XX, arma::trans(X) * y);
        arma::vec g_i(r, arma::fill::zeros);
        // Init g_incl_i using ridge regression to avoide zeros in the first iteration (invalid starting position)
        // TODO make this a parameter so that g_incl_i can be provided
        arma::vec g_incl_i = (Z.t() * Z + 10 * arma::eye(r, r)).i() * Z.t() * y;

        double s2_i = 1.0;

        // TODO make this a parameter so that w_i can be initialized with a given vector
        arma::Col<int> w_i(r, arma::fill::zeros);

        // TODO "Force include" vector as parameter -> geht als includevars in model selection

        arma::ivec z_cols_temp(z.n_cols);
        for (int i = 0; i < z.n_cols; ++i) {
            z_cols_temp(i) = i;
        }
        arma::ivec z_cols(z.n_cols * n);
        for (int i = 0; i < n; ++i) {
            z_cols.subvec(i * z.n_cols, (i + 1) * z.n_cols - 1) = z_cols_temp;
        }

        arma::vec w_1 = arma::conv_to<arma::vec>::from(z_cols) % w_i; // element-wise multiplication

        /* --------------------- Starting Values for HS-Plug-In --------------------- */
        arma::vec lamb_b(p, arma::fill::ones);
        double tau_b = 1.0;
        arma::vec nu_b(p, arma::fill::ones);
        double xi_b = 1.0;

        /* --------------------------- Prep Calculations ---------------------------- */
        double cN        = c0 + N / 2.0;
        arma::mat XX_inv = arma::inv(XX);

        // TODO this is never used. Comment it out for now
        // arma::vec Pr_g(n * (t - 2)); // NA in R can be represented by uninitialized values or a specific value if needed


        /* --------------------------------- Store ---------------------------------- */
        long Nstore = 0;
        if (num_burnin >= num_draws) {
            std::cerr << "Error: The number of burn-in exceeds number of draws" << std::endl;
            // Consider throwing an exception or handling this error appropriately in your program
            // TODO make function return early
            // return 0; // Or exit the function
        } else {
            Nstore = num_draws - num_burnin;
        }

        /* ---------------------------------- Main ---------------------------------- */
        arma::mat b_store(Nstore, p);
        arma::mat g_store(Nstore, r);
        arma::mat s2_store(Nstore, 1);

        arma::Mat<int> w_store(Nstore, r);

        // TODO never used. Commented out for now
        // arma::vec g_i_n0 = g_i;


        /* -------------------------------------------------------------------------- */
        /*                               Gibbs Sampler                                */
        /* -------------------------------------------------------------------------- */
        std::cout << "Starting Gibbs Sampler" << std::endl;


        for (int iter = 1 - num_burnin; iter <= Nstore; ++iter) {
            /* --------------------------- draw p(s2|a,b,g,y) --------------------------- */
            if (use_phiinit) {
                double CN = C0 + 0.5 * arma::sum(arma::square(y - X * b_i - Z * g_i));
                s2_i      = 1.0 / arma::randg(arma::distr_param(cN, 1.0 / CN)); // 1.0 / CN -> rate = 1/scale
            }

            /* --------------------------- draw p(b|a,g,s2,y) --------------------------- */
            if (b_prior == "hs") {
                // draw p(xi,nu|b,s2,y) ========.====
                arma::vec nu_b(lamb_b.n_elem);
                for (int j = 0; j < lamb_b.n_elem; ++j) {
                    nu_b(j) = 1.0 / arma::randg(arma::distr_param(1.0, 1.0 + 1.0 / lamb_b(j)));
                    if (nu_b(j) > 1e+8) nu_b(j) = 1e+8;
                    if (nu_b(j) < 1e-8) nu_b(j) = 1e-8;
                }
                xi_b = 1.0 / arma::randg(arma::distr_param(1.0, 1.0 + 1.0 / tau_b));
                if (xi_b > 1e+8) xi_b = 1e+8;
                if (xi_b < 1e-8) xi_b = 1e-8;

                // draw p(tau,lamb|xi,nu,b,s2,y).====
                arma::vec lamb_b_new(nu_b.n_elem);
                for (int j = 0; j < nu_b.n_elem; ++j) {
                    lamb_b_new(j) = 1.0 / arma::randg(
                                        arma::distr_param(
                                            1.0, 1.0 / nu_b(j) + std::pow(b_i(j), 2) / (2.0 * tau_b * s2_i)));
                    if (lamb_b_new(j) > 1e+8) lamb_b_new(j) = 1e+8;
                    if (lamb_b_new(j) < 1e-8) lamb_b_new(j) = 1e-8;
                }
                lamb_b = lamb_b_new;

                tau_b = 1.0 / arma::randg(arma::distr_param((p + 1.0) / 2.0,
                                                            1.0 / xi_b + arma::sum(arma::square(b_i) / lamb_b) / (
                                                                2.0 * s2_i)));
                if (tau_b > 1e+8) tau_b = 1e+8;
                if (tau_b < 1e-8) tau_b = 1e-8;

                arma::mat A     = XX + (1.0 / tau_b) * arma::diagmat(1.0 / lamb_b);
                arma::mat A_inv = arma::inv(A);
                arma::mat BN    = s2_i * A_inv;
                arma::vec bN    = A_inv * (arma::trans(X) * (y - Z * g_i));

                arma::mat BN_chol;
                bool chol_success = arma::chol(BN_chol, BN);
                if (chol_success) {
                    b_i = bN + BN_chol.t() * arma::randn(p);
                } else {
                    std::cerr <<
                            "Warning: Cholesky decomposition failed. Using multivariate normal from another method." <<
                            std::endl;
                    b_i = arma::mvnrnd(bN, BN);
                }
            } else if (b_prior == "g" || b_prior == "f") {
                arma::mat BN_inv = (1.0 / s2_i + 1.0 / lambda_b) * XX;
                arma::mat BN     = (s2_i * lambda_b / (s2_i + lambda_b)) * XX_inv;
                arma::vec bN     = (1.0 / (s2_i + lambda_b)) * (
                                   lambda_b * XX_inv * (arma::trans(X) * (y - Z * g_i)) + s2_i * b0);

                arma::mat BN_chol;
                bool chol_success = arma::chol(BN_chol, BN);
                if (chol_success) {
#ifdef DEV
                    // Standalone mode: Use C++11 <random>
                    arma::vec random_std_normals(p);
                    for (arma::uword i = 0; i < p; ++i) {
                        // Use the generator and distribution defined outside the loop
                        random_std_normals(i) = cpp_normal_dist(cpp_rng_generator);
                    }
                    b_i = bN + BN_chol.t() * random_std_normals;

#else
                    // R package mode: Use arma::randn (which will use R's RNG)
                    // No explicit seeding needed here; R controls it via set.seed()
                    b_i = bN + BN_chol.t() * arma::randn(p);
#endif
                } else {
                    std::cerr <<
                            "Warning: Cholesky decomposition failed (g/f prior). Using multivariate normal from another method."
                            << std::endl;
                    b_i = arma::mvnrnd(bN, BN);
                }
            }

            // arma::vec fixed_random = {1, -1, 2};
            // b_i                    = fixed_random;

            // ==================== draw p(w|a,b,s2,y) ====================
            arma::vec y_hat = y - X * b_i;

            /* ---------------------------- Model Selection ----------------------------- */

            int nn = 2;

            arma::vec thinit;
            InitType initpar_type;
            if (iter == 1 - num_burnin) {
                initpar_type = InitType::AUTO;
            } else {
                thinit = g_incl_i;
            }

            // TODO Make everything that is hardcoded a parameter.
            // TODO tau = step variance
            //  * priorDelta -> betabinom or bbprior.
            // priorDelta -> binary 1 -> modelbinomprior     2 -> modelbbprior
            //      -> wenn 1 -> prDelta.p -> prDeltap
            //      -> wenn 2 -> prDelta.a & .b -> vector in parprDeltap
            //      -> use defaults otherwise (0.5) (1,1)

            // TODO we need:
            // method
            // hesstype
            // optimmethod
            // optim maxit
            // B
            // knownhi
            // r
            // logscale
            // alpha
            // lambda
            // include vars (force include)


            // prior skew = TAU! ES GIBT NUR IN TAU
            // prDeltap = prConstrp & same for parprDeltap
            arma::Col<int> post_sample = model_selection_with_strategy(y_hat,
                                                                       Z,
                                                                       nn,
                                                                       1,
                                                                       nn - 1,
                                                                       w_i,
                                                                       false,
                                                                       false,
                                                                       true,
                                                                       s2_i,
                                                                       tau,
                                                                       0.348,
                                                                       0.5,
                                                                       thinit,
                                                                       initpar_type,
                                                                       new_par_method,
                                                                       new_par_hesstype,
                                                                       new_par_optim_method,
                                                                       new_par_optim_maxit,
                                                                       new_par_B,
                                                                       new_par_knownphi,
                                                                       new_par_r,
                                                                       new_par_alpha,
                                                                       new_par_lambda,
                                                                       strategy,
                                                                       n
            );

            w_i = post_sample;


            // ==================== draw p(g|w,a,b,s2,y) ====================
            if (arma::sum(w_i) > 0) {
                arma::vec ans = arma::vec().zeros(Z.n_cols);

                int nonZeroCount = arma::accu(w_i != 0);

                // Create a new matrix with the same number of rows as Z and columns equal to nonZeroCount
                arma::dmat filteredZ(Z.n_rows, nonZeroCount);

                // arma::vec filteredGi(nonZeroCount);

                int currentCol = 0;
                for (int i = 0; i < w_i.size(); ++i) {
                    if (w_i(i) != 0) {
                        // Copy the corresponding column from Z to filteredZ
                        filteredZ.col(currentCol) = Z.col(i);
                        // filteredGi(currentCol)    = g_i_n0(i);
                        ++currentCol;
                    }
                }

                // Filtered Z -> Nur die spalten wo w_i == 1

                bool use_thinit = false;
                if (iter > 1 - num_burnin) {
                    use_thinit = true;
                    // use_phiinit = true;
                }

                // TODO can be parallelized probably
                rnlpPost_lm(ans.memptr(),
                            1,
                            0,
                            1,
                            y_hat.memptr(),
                            filteredZ.memptr(),
                            y_hat.size(),
                            filteredZ.n_cols,
                            1,
                            tau_b,
                            c0,
                            C0,
                            1,
                            g_incl_i,
                            use_thinit,
                            s2_i,
                            use_phiinit);


                if (use_phiinit == false) {
                    s2_i = ans(nonZeroCount);
                }


                // TODO Document this!
                g_i = arma::vec().zeros(w_i.n_elem);

                int ans_pos = 0; // Position tracker for the ans vector
                for (int i = 0; i < w_i.n_elem; ++i) {
                    if (w_i(i) == 1) {
                        // When w_i is 1, use the next value from ans
                        g_i(i)      = ans[ans_pos];
                        g_incl_i(i) = ans[ans_pos];
                        ans_pos++; // Only increment ans_pos when we use a value
                    }
                    // When w_i is 0, g_i stays 0 (already initialized that way)
                }

                // std::cout << "w_i:\n" << w_i << std::endl;
                // std::cout << "g_i:\n" << g_i << std::endl;
            } else {
                g_i.zeros(r);
            }

            // TODO is this a correct nan check?
            if (!arma::find_nan(g_i).is_empty()) {
                std::cerr << "NaNs Produced in g_i at iteration " << iter << std::endl;
            }

            // ==================== store and adjust tau ====================
            if (iter > 0) {
                // R's loop goes from 1-Nburn, so storing starts after burn-in
                b_store.row(iter - 1) = b_i.t();
                g_store.row(iter - 1) = g_i.t();
                w_store.row(iter - 1) = w_i.t();
                s2_store(iter - 1, 0) = s2_i;
            } else {
                // Adjust tau during burn-in (placeholder logic)
                // if (arma::sum(g_i) < n * t * 0.05) {
                //     tau = tau * 5.0 - arma::sum(g_i);
                // }
            }
        }


        /* -------------------------------------------------------------------------- */
        /*                        Calculate and Print w_store Means                      */
        /* -------------------------------------------------------------------------- */
        std::cout << "Calculating column means of w_store..." << std::endl;

        // Calculate column means of w_store
        arma::rowvec w_store_means(r);
        for (int j = 0; j < r; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Nstore; ++i) {
                sum += w_store(i, j);
            }
            w_store_means(j) = sum / Nstore;
        }

        // Calculate column means of w_store
        arma::rowvec b_store_means(p);
        for (int j = 0; j < p; ++j) {
            double sum = 0.0;
            for (int i = 0; i < Nstore; ++i) {
                sum += b_store(i, j);
            }
            b_store_means(j) = sum / Nstore;
        }

        arma::rowvec s2_store_means(1);
        double sum = 0.0;
        for (int i = 0; i < Nstore; ++i) {
            sum += s2_store(i, 0);
        }
        s2_store_means(0) = sum / Nstore;

        return BisamResult{
            b_store, g_store, s2_store, w_store, w_store_means, b_store_means, s2_store_means
        };
    }
}

