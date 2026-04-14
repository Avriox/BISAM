#ifndef SIMPLE_RESULT_STORAGE_H
#define SIMPLE_RESULT_STORAGE_H

#include <string>
#include <vector>
#include <fstream>
#include <filesystem>
#include <chrono>
#include <iomanip>
#include <sstream>
#include <RcppArmadillo.h>
#include "../include/biasm_model.h"
#include "../r-testing-files/simulation_data/simulation_datasets.h"
#include "validation.h"

namespace bisam {
    struct TimingResult {
        std::string experiment_name;
        std::string run_name; // NEW: identifies this specific run/version
        std::string dataset_name;
        int run_number;
        double execution_time_ms;
        int n;
        int t;
        int nx;
        std::string timestamp;

        TimingResult(const std::string &exp_name, const std::string &r_name,
                     const std::string &ds_name, int run_num, double time_ms,
                     int n_val, int t_val, int nx_val)
            : experiment_name(exp_name), run_name(r_name), dataset_name(ds_name),
              run_number(run_num), execution_time_ms(time_ms), n(n_val), t(t_val), nx(nx_val) {
            // Add timestamp
            auto now    = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
            timestamp = ss.str();
        }
    };

    struct EstimationResult {
        std::string experiment_name;
        std::string run_name; // NEW: identifies this specific run/version
        std::string dataset_name;
        int run_number;
        std::string parameter_type; // "beta", "constant", "break_indicator", "sigma2"
        int parameter_index;
        double true_value;
        double estimated_value;
        double error;
        double probability; // For break indicators
        std::string timestamp;

        EstimationResult(const std::string &exp_name, const std::string &r_name,
                         const std::string &ds_name, int run_num, const std::string &param_type,
                         int param_idx, double true_val, double est_val, double err, double prob = 0.0)
            : experiment_name(exp_name), run_name(r_name), dataset_name(ds_name),
              run_number(run_num), parameter_type(param_type), parameter_index(param_idx),
              true_value(true_val), estimated_value(est_val), error(err), probability(prob) {
            auto now    = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
            timestamp = ss.str();
        }
    };

    struct ValidationSummary {
        std::string experiment_name;
        std::string run_name; // NEW: identifies this specific run/version
        std::string dataset_name;
        int run_number;
        double beta_mse;
        double beta_max_error;
        bool break_detected_correctly;
        int true_break_position;
        int detected_break_position;
        double detected_break_probability;
        double constant_error;
        double detection_threshold;
        std::string timestamp;

        ValidationSummary(const std::string &exp_name, const std::string &r_name,
                          const std::string &ds_name, int run_num, const ValidationResults &val)
            : experiment_name(exp_name), run_name(r_name), dataset_name(ds_name),
              run_number(run_num), beta_mse(val.beta_mse), beta_max_error(val.beta_max_error),
              break_detected_correctly(val.break_detected_correctly),
              true_break_position(val.true_break_position),
              detected_break_position(val.detected_break_position),
              detected_break_probability(val.detected_break_probability),
              constant_error(val.constant_error),
              detection_threshold(val.detection_threshold) {
            auto now    = std::chrono::system_clock::now();
            auto time_t = std::chrono::system_clock::to_time_t(now);
            std::stringstream ss;
            ss << std::put_time(std::localtime(&time_t), "%Y-%m-%d %H:%M:%S");
            timestamp = ss.str();
        }
    };

    class ResultsStorage {
    private:
        std::string experiment_name_;
        std::string base_directory_;
        bool append_mode_;
        std::vector<TimingResult> timing_results_;
        std::vector<EstimationResult> estimation_results_;
        std::vector<ValidationSummary> validation_summaries_;

        static arma::vec compute_beta_means(const arma::mat &beta_samples) {
            return arma::mean(beta_samples, 0).t();
        }

        void ensure_directory_exists(const std::string &path) {
            if (!std::filesystem::exists(path)) {
                std::filesystem::create_directories(path);
            }
        }

        void delete_file_if_exists(const std::string &filepath) {
            if (std::filesystem::exists(filepath)) {
                std::filesystem::remove(filepath);
            }
        }

        bool file_exists_and_has_content(const std::string &filepath) {
            if (!std::filesystem::exists(filepath)) return false;
            std::ifstream file(filepath);
            return file.good() && file.peek() != std::ifstream::traits_type::eof();
        }

    public:
        ResultsStorage(const std::string &experiment_name, bool append_mode = false,
                       const std::string &base_dir                          = "simulations")
            : experiment_name_(experiment_name), base_directory_(base_dir), append_mode_(append_mode) {
            ensure_directory_exists(base_directory_);

            if (!append_mode) {
                // Delete existing files if not in append mode
                std::vector<std::string> file_patterns = {
                    experiment_name_ + "_timing.csv",
                    experiment_name_ + "_estimates.csv",
                    experiment_name_ + "_validation.csv"
                };

                for (const auto &pattern: file_patterns) {
                    delete_file_if_exists(base_directory_ + "/" + pattern);
                }
            }
        }

        void add_timing_result(const std::string &dataset_name, const std::string &run_name,
                               int run_number, double execution_time_ms, int n, int t, int nx) {
            timing_results_.emplace_back(experiment_name_, run_name, dataset_name, run_number,
                                         execution_time_ms, n, t, nx);
        }

        void add_estimation_results(const std::string &dataset_name, const std::string &run_name,
                                    int run_number, const bisam::BisamResult &result, const DatasetInfo &dataset) {
            // Store beta estimates
            if (result.beta_samples.n_cols > 0 && dataset.true_beta.n_elem > 0) {
                arma::vec estimated_betas = compute_beta_means(result.beta_samples);

                // Only process matching dimensions
                int n_betas = std::min(static_cast<int>(estimated_betas.n_elem),
                                       static_cast<int>(dataset.true_beta.n_elem));

                for (int i = 0; i < n_betas; i++) {
                    double error = std::abs(estimated_betas(i) - dataset.true_beta(i));
                    estimation_results_.emplace_back(experiment_name_, run_name, dataset_name, run_number,
                                                     "beta", i, dataset.true_beta(i),
                                                     estimated_betas(i), error);
                }
            }

            // Store constant estimate if present
            if (dataset.has_const && result.beta_samples.n_cols > dataset.true_beta.n_elem) {
                arma::vec estimated_betas = compute_beta_means(result.beta_samples);
                double estimated_const    = estimated_betas(estimated_betas.n_elem - 1);
                double error              = std::abs(estimated_const - dataset.true_const);
                estimation_results_.emplace_back(experiment_name_, run_name, dataset_name, run_number,
                                                 "constant", 0, dataset.true_const,
                                                 estimated_const, error);
            }

            // Store break indicators
            // RUNTIME DIMENSIONS: with z = z_lower_tri.cols(2, t-2), r = n*(t-3)
            int R = static_cast<int>(result.indicator_means.n_elem);

            // Precompute truth flags at model Z indices (0-based)
            std::vector<char> is_true_break(R, 0);
            if (dataset.true_breaks.n_rows > 0) {
                for (int j = 0; j < dataset.true_breaks.n_rows; ++j) {
                    // Column 1 stores the 0-based Z column index as double; round safely before casting
                    long long zidx_ll = static_cast<long long>(std::llround(dataset.true_breaks(j, 1)));
                    if (zidx_ll >= 0 && zidx_ll < R) {
                        is_true_break[static_cast<int>(zidx_ll)] = 1;
                    }
                }
            }

            for (int i = 0; i < R; i++) {
                double true_indicator = is_true_break[i] ? 1.0 : 0.0;
                double est            = result.indicator_means(i);
                double error          = std::abs(est - true_indicator);

                estimation_results_.emplace_back(experiment_name_, run_name, dataset_name, run_number,
                                                 "break_indicator", i, true_indicator,
                                                 est, error, est);
            }

            // Store sigma2 estimate
            if (result.sigma2_means.n_elem > 0) {
                estimation_results_.emplace_back(experiment_name_, run_name, dataset_name, run_number,
                                                 "sigma2", 0, dataset.error_sd * dataset.error_sd,
                                                 result.sigma2_means(0),
                                                 std::abs(
                                                     result.sigma2_means(0) - dataset.error_sd * dataset.error_sd));
            }
        }

        void add_validation_summary(const std::string &dataset_name, const std::string &run_name,
                                    int run_number, const ValidationResults &validation) {
            validation_summaries_.emplace_back(experiment_name_, run_name, dataset_name, run_number, validation);
        }

        // Replace the save_results() method in your ResultsStorage class with this debug version:

        void save_results() {
            std::cout << "=== SAVE_RESULTS DEBUG START ===" << std::endl;

            // Get absolute path
            std::filesystem::path abs_path = std::filesystem::absolute(base_directory_);
            std::cout << "Base directory: " << abs_path << std::endl;
            std::cout << "Directory exists: " << (std::filesystem::exists(abs_path) ? "YES" : "NO") << std::endl;

            // Check data sizes
            std::cout << "Data to save:" << std::endl;
            std::cout << "  - Timing results: " << timing_results_.size() << " entries" << std::endl;
            std::cout << "  - Estimation results: " << estimation_results_.size() << " entries" << std::endl;
            std::cout << "  - Validation summaries: " << validation_summaries_.size() << " entries" << std::endl;

            if (timing_results_.empty() && estimation_results_.empty() && validation_summaries_.empty()) {
                std::cout << "WARNING: No data to save!" << std::endl;
            }

            try {
                std::cout << "Attempting to save timing results..." << std::endl;
                save_timing_results();
                std::cout << "✓ Timing results saved successfully" << std::endl;
            } catch (const std::exception &e) {
                std::cout << "✗ ERROR saving timing results: " << e.what() << std::endl;
            }

            try {
                std::cout << "Attempting to save estimation results..." << std::endl;
                save_estimation_results();
                std::cout << "✓ Estimation results saved successfully" << std::endl;
            } catch (const std::exception &e) {
                std::cout << "✗ ERROR saving estimation results: " << e.what() << std::endl;
            }

            try {
                std::cout << "Attempting to save validation summaries..." << std::endl;
                save_validation_summaries();
                std::cout << "✓ Validation summaries saved successfully" << std::endl;
            } catch (const std::exception &e) {
                std::cout << "✗ ERROR saving validation summaries: " << e.what() << std::endl;
            }

            std::cout << "=== SAVE_RESULTS DEBUG END ===" << std::endl;
        }

        // Also replace save_timing_results() with this debug version:
        void save_timing_results() {
            std::string filename               = base_directory_ + "/" + experiment_name_ + "_timing.csv";
            std::filesystem::path abs_filename = std::filesystem::absolute(filename);

            std::cout << "  Timing file: " << abs_filename << std::endl;
            std::cout << "  Data entries: " << timing_results_.size() << std::endl;

            if (timing_results_.empty()) {
                std::cout << "  WARNING: No timing data to save!" << std::endl;
                return;
            }

            bool should_append = append_mode_ && file_exists_and_has_content(filename);
            std::cout << "  Append mode: " << (should_append ? "YES" : "NO") << std::endl;

            std::ofstream file;
            if (should_append) {
                file.open(filename, std::ios::app);
                std::cout << "  Opened in append mode" << std::endl;
            } else {
                file.open(filename);
                std::cout << "  Opened in write mode" << std::endl;
                // Write header for new file
                file << "experiment_name,run_name,dataset_name,run_number,execution_time_ms,n,t,nx,timestamp\n";
                std::cout << "  Header written" << std::endl;
            }

            if (!file.is_open()) {
                throw std::runtime_error("Could not open timing results file: " + filename);
            }

            std::cout << "  File opened successfully, writing " << timing_results_.size() << " entries..." <<
                    std::endl;

            // Write data
            int count = 0;
            for (const auto &result: timing_results_) {
                file << result.experiment_name << ","
                        << result.run_name << ","
                        << result.dataset_name << ","
                        << result.run_number << ","
                        << std::fixed << std::setprecision(3) << result.execution_time_ms << ","
                        << result.n << ","
                        << result.t << ","
                        << result.nx << ","
                        << result.timestamp << "\n";
                count++;
            }

            std::cout << "  Wrote " << count << " entries" << std::endl;

            file.close();
            std::cout << "  File closed" << std::endl;

            // Verify file was created
            if (std::filesystem::exists(abs_filename)) {
                auto file_size = std::filesystem::file_size(abs_filename);
                std::cout << "  ✓ File created successfully, size: " << file_size << " bytes" << std::endl;
            } else {
                std::cout << "  ✗ File was not created!" << std::endl;
            }
        }

        // NEW: Clear memory after saving to prevent accumulation
        void clear_memory() {
            timing_results_.clear();
            estimation_results_.clear();
            validation_summaries_.clear();

            // Optional: Shrink to fit to release memory back to OS
            timing_results_.shrink_to_fit();
            estimation_results_.shrink_to_fit();
            validation_summaries_.shrink_to_fit();
        }

    private
    :
        // void save_timing_results() {
        //     std::string filename = base_directory_ + "/" + experiment_name_ + "_timing.csv";
        //     bool should_append = append_mode_ && file_exists_and_has_content(filename);
        //
        //     std::ofstream file;
        //     if (should_append) {
        //         file.open(filename, std::ios::app);
        //     } else {
        //         file.open(filename);
        //         // Write header for new file
        //         file << "experiment_name,run_name,dataset_name,run_number,execution_time_ms,n,t,nx,timestamp\n";
        //     }
        //
        //     if (!file.is_open()) {
        //         throw std::runtime_error("Could not open timing results file: " + filename);
        //     }
        //
        //     // Write data
        //     for (const auto& result : timing_results_) {
        //         file << result.experiment_name << ","
        //              << result.run_name << ","
        //              << result.dataset_name << ","
        //              << result.run_number << ","
        //              << std::fixed << std::setprecision(3) << result.execution_time_ms << ","
        //              << result.n << ","
        //              << result.t << ","
        //              << result.nx << ","
        //              << result.timestamp << "\n";
        //     }
        //
        //     file.close();
        // }

        void save_estimation_results() {
            std::string filename = base_directory_ + "/" + experiment_name_ + "_estimates.csv";
            bool should_append   = append_mode_ && file_exists_and_has_content(filename);

            std::ofstream file;
            if (should_append) {
                file.open(filename, std::ios::app);
            } else {
                file.open(filename);
                file << "experiment_name,run_name,dataset_name,run_number,parameter_type,"
                        << "parameter_index,true_value,estimated_value,error,probability,timestamp\n";
            }

            if (!file.is_open()) {
                throw std::runtime_error("Could not open estimation results file: " + filename);
            }

            // Write data
            for (const auto &result: estimation_results_) {
                file << result.experiment_name << ","
                        << result.run_name << ","
                        << result.dataset_name << ","
                        << result.run_number << ","
                        << result.parameter_type << ","
                        << result.parameter_index << ","
                        << std::fixed << std::setprecision(8) << result.true_value << ","
                        << std::setprecision(8) << result.estimated_value << ","
                        << std::setprecision(8) << result.error << ","
                        << std::setprecision(8) << result.probability << ","
                        << result.timestamp << "\n";
            }

            file.close();
        }

        void save_validation_summaries() {
            std::string filename = base_directory_ + "/" + experiment_name_ + "_validation.csv";
            bool should_append   = append_mode_ && file_exists_and_has_content(filename);

            std::ofstream file;
            if (should_append) {
                file.open(filename, std::ios::app);
            } else {
                file.open(filename);
                file << "experiment_name,run_name,dataset_name,run_number,beta_mse,beta_max_error,"
                        << "break_detected_correctly,true_break_position,detected_break_position,"
                        << "detected_break_probability,constant_error,detection_threshold,timestamp\n";
            }

            if (!file.is_open()) {
                throw std::runtime_error("Could not open validation results file: " + filename);
            }

            // Write data
            for (const auto &summary: validation_summaries_) {
                file << summary.experiment_name << ","
                        << summary.run_name << ","
                        << summary.dataset_name << ","
                        << summary.run_number << ","
                        << std::fixed << std::setprecision(8) << summary.beta_mse << ","
                        << std::setprecision(8) << summary.beta_max_error << ","
                        << (summary.break_detected_correctly ? "TRUE" : "FALSE") << ","
                        << summary.true_break_position << ","
                        << summary.detected_break_position << ","
                        << std::setprecision(8) << summary.detected_break_probability << ","
                        << std::setprecision(8) << summary.constant_error << ","
                        << std::setprecision(8) << summary.detection_threshold << ","
                        << summary.timestamp << "\n";
            }

            file.close();
        }
    };
} // namespace bisam

#endif // SIMPLE_RESULT_STORAGE_H
