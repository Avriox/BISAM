//
// Created by jakob on 4/15/25.
//

#ifndef FAST_POLYNOMIAL_SOLVER_H
#define FAST_POLYNOMIAL_SOLVER_H

#include <algorithm>
#include <functional>
#include <cfloat>
#include <unordered_map>
#include <cstdint>
#include <cstring>
#include <cmath>

// The original code used _DBL_RADIX which is not defined
# define _DBL_RADIX FLT_RADIX

using namespace std;

// ===== Common constants for both solvers =====
constexpr int POLYNOMIAL_DEGREE = 4;                     // We always have degree 4 polynomials
constexpr int MAX_COEFF         = POLYNOMIAL_DEGREE + 1; // Number of coefficients in a degree 4 polynomial
constexpr int MAX_ROOTS         = POLYNOMIAL_DEGREE;     // Maximum number of roots

// ===== Newton solver constants =====
constexpr int MAX_ITER = 50;

// ===== Aberth-Ehrlich solver constants =====
constexpr int MAX_ABERTH_ITER          = 50;
constexpr double OSCILLATION_THRESHOLD = 1e-10; // Threshold to detect oscillation
constexpr double GOOD_ENOUGH_THRESHOLD = 1e-15; // Early termination threshold
constexpr int MAX_OSCILLATION_COUNT    = 4;     // Maximum oscillation count before early termination

// ===== Global cache and helper structs =====

// Structure to store roots with a count of valid roots
struct RootStorage {
    double real[MAX_ROOTS]; // Real parts of the roots
    double imag[MAX_ROOTS]; // Imaginary parts of the roots
    int count;              // Number of valid roots
};

// Thread-local cache for storing polynomial roots based on coefficient hash
inline thread_local std::unordered_map<uint32_t, RootStorage> g_rootCache;

// Structure representing a complex number without using std::complex
struct Complex {
    double re; // Real part
    double im; // Imaginary part

    // Fast norm (squared magnitude)
    inline double norm() const {
        return re * re + im * im;
    }

    // Fast absolute value (magnitude)
    inline double abs() const {
        return sqrt(norm());
    }
};

// Structure for polynomial evaluation
struct Eval {
    double z_re;  // Real part of z
    double z_im;  // Imaginary part of z
    double pz_re; // Real part of P(z)
    double pz_im; // Imaginary part of P(z)
    double apz;   // |P(z)| - absolute value
};

// ===== Complex arithmetic operations =====

// Complex number addition
inline void complex_add(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    result_re = a_re + b_re;
    result_im = a_im + b_im;
}

// Complex number subtraction
inline void complex_sub(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    result_re = a_re - b_re;
    result_im = a_im - b_im;
}

// Complex number multiplication
inline void complex_mul(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    double tmp_re = a_re * b_re - a_im * b_im;
    double tmp_im = a_re * b_im + a_im * b_re;
    result_re     = tmp_re;
    result_im     = tmp_im;
}

// Complex number division
inline void complex_div(double &result_re, double &result_im,
                        double a_re, double a_im,
                        double b_re, double b_im) {
    double denom  = b_re * b_re + b_im * b_im;
    double tmp_re = (a_re * b_re + a_im * b_im) / denom;
    double tmp_im = (a_im * b_re - a_re * b_im) / denom;
    result_re     = tmp_re;
    result_im     = tmp_im;
}

// Complex number absolute value (magnitude)
inline double complex_abs(double re, double im) {
    return sqrt(re * re + im * im);
}

// Complex number conjugate
inline void complex_conj(double &result_re, double &result_im,
                         double a_re, double a_im) {
    result_re = a_re;
    result_im = -a_im;
}

// Scale complex number by a real value
inline void complex_scale(double &result_re, double &result_im,
                          double a_re, double a_im,
                          double scale) {
    result_re = a_re * scale;
    result_im = a_im * scale;
}

// ===== Common helper function declarations =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const double coefficients[MAX_COEFF]);

// Root cache management functions
static void ClearRootCache();

static void PrintRootCacheStats();

static void ReserveRootCache(size_t expectedSize);

// ===== Newton solver helper function declarations =====

// Compute suitable starting point for root finding
static double ComputeStartPoint(const double a[MAX_COEFF], int n);

// Evaluate polynomial with Horner's method
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, double z_re, double z_im);

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, double z_re, double z_im);

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, double x);

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, double z_re, double z_im);

// Solve quadratic polynomial
static int SolveQuadratic(const double a[MAX_COEFF], int n,
                          double real_roots[MAX_ROOTS], double imag_roots[MAX_ROOTS],
                          int rootIndex, double new_real[MAX_ROOTS], double new_imag[MAX_ROOTS],
                          int newRootIndex);

// Main function to find polynomial roots using Newton's method
int PolynomialRootsNewton(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr);

// ===== Aberth-Ehrlich solver helper function declarations =====

// Helper function to eliminate zero roots (roots that are exactly zero)
static int ZeroRoots(int n, double a_re[], double a_im[], double res_re[], double res_im[]);

// Helper function to evaluate a polynomial at a complex point using Horner's method
static Eval Horner(int n, const double a_re[], const double a_im[], double z_re, double z_im);

// Helper function to calculate starting points for the algorithm
static void StartPoints(int n, const double apolyr[], double Z_re[], double Z_im[]);

// Helper function to solve quadratic and linear equations directly
static void Quadratic(int n, double a_re[], double a_im[], double res_re[], double res_im[]);

// Main implementation of the Aberth-Ehrlich method
static void AberthEhrlich(int n, const double coeff_re[], const double coeff_im[],
                          double res_re[], double res_im[],
                          const double initialGuesses_re[] = nullptr,
                          const double initialGuesses_im[] = nullptr,
                          int numGuesses                   = 0);

// Main function to find polynomial roots using Aberth-Ehrlich method
int PolynomialRootsAberth(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr);


#endif //FAST_POLYNOMIAL_SOLVER_H
