#ifndef VALIDATION_H
#define VALIDATION_H

#include <RcppArmadillo.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include "../include/biasm_model.h"
#include "../r-testing-files/simulation_data/simulation_datasets.h"

struct ValidationResults {
    double beta_mse;
    double beta_max_error;
    bool break_detected_correctly;
    int true_break_position;
    int detected_break_position;
    double detected_break_probability;
    double break_magnitude_error;
    double detection_threshold;
    bool has_false_positive;
    double constant_error;
};

class BisamValidator {
public:
    static ValidationResults validate_results(const bisam::BisamResult &result,
                                              const DatasetInfo &dataset,
                                              double detection_threshold = 0.5,
                                              bool print_output          = true) {
        ValidationResults val;
        val.detection_threshold        = detection_threshold;
        val.beta_mse                   = 0.0;
        val.beta_max_error             = 0.0;
        val.constant_error             = 0.0;
        val.true_break_position        = -1;
        val.detected_break_position    = -1;
        val.detected_break_probability = 0.0;
        val.break_detected_correctly   = false;
        val.has_false_positive         = false;
        val.break_magnitude_error      = 0.0;

        validate_betas(val, result, dataset, print_output);
        validate_constant(val, result, dataset, print_output);
        validate_breaks(val, result, dataset, print_output);

        return val;
    }

    static void print_dataset_info(const DatasetInfo &dataset) {
        std::cout << "\n" << std::string(70, '=') << std::endl;
        std::cout << "DATASET: " << dataset.name
                << " (n=" << dataset.n << ", t=" << dataset.t << ", nx=" << dataset.nx << ")" << std::endl;
        std::cout << std::string(70, '=') << std::endl;

        std::cout << "Simulation parameters:" << std::endl;
        std::cout << "  Error SD: " << dataset.error_sd << std::endl;
        std::cout << "  Step mean: " << dataset.step_mean << std::endl;
        std::cout << "  Has constant: " << (dataset.has_const ? "Yes" : "No") << std::endl;

        if (dataset.has_const) {
            std::cout << "  True constant: " << dataset.true_const << std::endl;
        }

        std::cout << "\nData info:" << std::endl;
        std::cout << "  Dimensions: " << dataset.data.n_rows << " x " << dataset.data.n_cols << std::endl;
        std::cout << "  True betas: ";
        if (dataset.true_beta.n_elem > 0) {
            dataset.true_beta.t().raw_print(std::cout);
        } else {
            std::cout << "None" << std::endl;
        }

        if (dataset.true_breaks.n_rows > 0) {
            std::cout << "  True breaks: " << dataset.true_breaks.n_rows << " break(s)" << std::endl;
        } else {
            std::cout << "  True breaks: None" << std::endl;
        }
    }

    static void print_validation_summary(const ValidationResults &val, const std::string &dataset_name) {
        std::cout << "\n" << std::string(70, '-') << std::endl;
        std::cout << "VALIDATION SUMMARY for " << dataset_name << std::endl;
        std::cout << std::string(70, '-') << std::endl;

        // Beta validation
        if (val.beta_max_error > 0) {
            std::cout << "Beta estimation: "
                    << (val.beta_max_error < 0.5 ? "GOOD" : "NEEDS ATTENTION")
                    << " (MSE: " << std::scientific << std::setprecision(3) << val.beta_mse
                    << ", Max error: " << std::fixed << std::setprecision(4) << val.beta_max_error << ")" << std::endl;
        }

        // Constant validation
        if (val.constant_error > 0) {
            std::cout << "Constant estimation: "
                    << (val.constant_error < 0.5 ? "GOOD" : "NEEDS ATTENTION")
                    << " (Error: " << val.constant_error << ")" << std::endl;
        }

        // Break validation
        std::cout << "Break detection: ";
        if (val.true_break_position >= 0) {
            std::cout << (val.break_detected_correctly ? "CORRECT" : "INCORRECT")
                    << " (Detected prob: " << val.detected_break_probability << ")" << std::endl;
            std::cout << "  Position mapping: true=" << val.true_break_position
                    << " -> detected=" << val.detected_break_position << std::endl;
        } else {
            std::cout << (val.break_detected_correctly ? "CORRECT (no false positives)" : "FALSE POSITIVE")
                    << " (Max prob: " << val.detected_break_probability << ")" << std::endl;
        }

        std::cout << std::string(70, '-') << std::endl;
    }

    static void print_detailed_results(const bisam::BisamResult &result, const DatasetInfo &dataset) {
        std::cout << "\n=== DETAILED MCMC RESULTS ===" << std::endl;

        std::cout << "MCMC Summary:" << std::endl;
        std::cout << "  Beta samples: " << result.beta_samples.n_rows << " samples x "
                << result.beta_samples.n_cols << " parameters" << std::endl;
        std::cout << "  Sigma2 mean: " << std::fixed << std::setprecision(4) << result.sigma2_means(0) << std::endl;

        // Show first few indicator probabilities
        std::cout << "\nFirst 15 indicator probabilities:" << std::endl;
        int n_to_show = std::min(15, static_cast<int>(result.indicator_means.n_elem));
        for (int i = 0; i < n_to_show; i++) {
            std::cout << std::setw(8) << std::fixed << std::setprecision(4) << result.indicator_means(i);
            if ((i + 1) % 5 == 0) std::cout << std::endl;
        }
        if (n_to_show % 5 != 0) std::cout << std::endl;

        // Show any high-probability indicators
        std::cout << "\nHigh-probability indicators (> 0.1):" << std::endl;
        bool found_high = false;
        for (int i = 0; i < result.indicator_means.n_elem; i++) {
            if (result.indicator_means(i) > 0.1) {
                std::cout << "  Position " << i << ": " << result.indicator_means(i) << std::endl;
                found_high = true;
            }
        }
        if (!found_high) {
            std::cout << "  None found" << std::endl;
        }
    }

private:
    static arma::vec compute_beta_means(const arma::mat &beta_samples) {
        // beta_samples is (n_samples x n_betas), we want column means
        return arma::mean(beta_samples, 0).t(); // transpose to get column vector
    }

    static void validate_betas(ValidationResults &val, const bisam::BisamResult &result, const DatasetInfo &dataset,
                               bool print_output) {
        if (dataset.true_beta.n_elem > 0 && result.beta_samples.n_cols > 0) {
            arma::vec estimated_betas = compute_beta_means(result.beta_samples);

            // Ensure dimensions match
            if (estimated_betas.n_elem == dataset.true_beta.n_elem) {
                arma::vec beta_errors = arma::abs(estimated_betas - dataset.true_beta);
                val.beta_mse          = arma::mean(arma::square(beta_errors));
                val.beta_max_error    = arma::max(beta_errors);

                if (print_output) {
                    std::cout << "\n=== BETA VALIDATION ===" << std::endl;
                    std::cout << "True betas:      ";
                    dataset.true_beta.t().raw_print(std::cout);
                    std::cout << "Estimated betas: ";
                    estimated_betas.t().raw_print(std::cout);
                    std::cout << "Errors:          ";
                    beta_errors.t().raw_print(std::cout);
                    std::cout << "Beta MSE: " << std::fixed << std::setprecision(6) << val.beta_mse << std::endl;
                    std::cout << "Beta Max Error: " << val.beta_max_error << std::endl;
                }
            } else {
                if (print_output) {
                    std::cout << "\n=== BETA VALIDATION ERROR ===" << std::endl;
                    std::cout << "Dimension mismatch - True: " << dataset.true_beta.n_elem
                            << ", Estimated: " << estimated_betas.n_elem << std::endl;
                }
            }
        }
    }

    static void validate_constant(ValidationResults &val, const bisam::BisamResult &result,
                                  const DatasetInfo &dataset, bool print_output) {
        if (dataset.has_const) {
            // Assuming constant is the last column in beta_samples if present
            if (result.beta_samples.n_cols > dataset.true_beta.n_elem) {
                arma::vec estimated_betas = compute_beta_means(result.beta_samples);
                double estimated_const    = estimated_betas(estimated_betas.n_elem - 1);
                val.constant_error        = std::abs(estimated_const - dataset.true_const);

                if (print_output) {
                    std::cout << "\n=== CONSTANT VALIDATION ===" << std::endl;
                    std::cout << "True constant: " << dataset.true_const << std::endl;
                    std::cout << "Estimated constant: " << estimated_const << std::endl;
                    std::cout << "Constant error: " << val.constant_error << std::endl;
                }
            }
        }
    }

    static void validate_breaks(ValidationResults &val, const bisam::BisamResult &result, const DatasetInfo &dataset,
                                bool print_output) {
        // Highest indicator probability and its index
        if (result.indicator_means.n_elem == 0) return;
        val.detected_break_probability = arma::max(result.indicator_means);
        val.detected_break_position    = arma::index_max(result.indicator_means);

        if (print_output) {
            std::cout << "\n=== BREAK VALIDATION ===" << std::endl;
        }

        // Collect true Z indices directly from dataset.true_breaks
        std::vector<int> true_z_positions;
        if (dataset.true_breaks.n_rows > 0) {
            for (int j = 0; j < dataset.true_breaks.n_rows; ++j) {
                long long zidx_ll = static_cast<long long>(std::llround(dataset.true_breaks(j, 1)));
                if (zidx_ll >= 0 && zidx_ll < static_cast<long long>(result.indicator_means.n_elem)) {
                    true_z_positions.push_back(static_cast<int>(zidx_ll));
                }
            }
        }

        // For summary fields, record first true position (if any)
        val.true_break_position = true_z_positions.empty() ? -1 : true_z_positions.front();

        if (!true_z_positions.empty()) {
            if (print_output) {
                std::cout << "True break Z indices (0-based): ";
                for (size_t k = 0; k < true_z_positions.size(); ++k) {
                    std::cout << true_z_positions[k] << (k + 1 < true_z_positions.size() ? ", " : "");
                }
                std::cout << std::endl;
                std::cout << "Detected break: position " << val.detected_break_position
                        << " (prob: " << std::fixed << std::setprecision(4)
                        << val.detected_break_probability << ")" << std::endl;
            }

            bool detected_matches_truth = false;
            for (int zidx: true_z_positions) {
                if (val.detected_break_position == zidx) {
                    detected_matches_truth = true;
                    break;
                }
            }

            val.break_detected_correctly =
                    detected_matches_truth && (val.detected_break_probability > val.detection_threshold);

            // Optional placeholder for magnitude error
            val.break_magnitude_error = 0.0;
        } else {
            // No true breaks
            if (print_output) {
                std::cout << "No true breaks in dataset" << std::endl;
                std::cout << "Highest indicator probability: " << val.detected_break_probability
                        << " at position " << val.detected_break_position << std::endl;
            }
            val.has_false_positive       = (val.detected_break_probability > val.detection_threshold);
            val.break_detected_correctly = !val.has_false_positive;
        }

        // Show indicator probabilities around the detected break
        if (print_output) {
            std::cout << "\nIndicator probabilities around detected position:" << std::endl;
            int start = std::max(0, val.detected_break_position - 5);
            int end   = std::min(static_cast<int>(result.indicator_means.n_elem),
                               val.detected_break_position + 6);

            for (int i = start; i < end; i++) {
                std::cout << "  Position " << std::setw(2) << i << ": "
                        << std::fixed << std::setprecision(4) << result.indicator_means(i);
                if (i == val.detected_break_position) std::cout << " <-- DETECTED";
                for (int zidx: true_z_positions) {
                    if (i == zidx) {
                        std::cout << " <-- TRUE";
                        break;
                    }
                }
                std::cout << std::endl;
            }
        }
    }
};

#endif // VALIDATION_H
