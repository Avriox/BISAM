//
// Created by jakob on 4/15/25.
//


#include "fast_polynomial_solver.h"

#include <iostream>
#include <ostream>

// ===== Common helper functions =====

// Ultra-fast hashing function specialized for our specific polynomial structure
static uint32_t HashPolynomial(const double coefficients[MAX_COEFF]) {
    // Get the integer parts of coefficients at positions 4, 1, and 0
    int c4_int = static_cast<int>(coefficients[4]);
    int c1_int = static_cast<int>(coefficients[1]);
    int c0_int = static_cast<int>(coefficients[0]);

    // Convert to positive values for consistent hashing
    uint32_t c4 = static_cast<uint32_t>(c4_int < 0 ? -c4_int : c4_int);
    uint32_t c1 = static_cast<uint32_t>(c1_int < 0 ? -c1_int : c1_int);
    uint32_t c0 = static_cast<uint32_t>(c0_int < 0 ? -c0_int : c0_int);

    // Create hash by bit-packing (very fast)
    // We reserve 10 bits for each value (allows values up to 1023)
    // Add sign bits in highest positions
    uint32_t hash = (c0 << 20) | (c1 << 10) | c4;
    // Include sign information in the two highest bits
    if (c0_int < 0) hash |= (1U << 30);
    if (c1_int < 0) hash |= (1U << 31);

    return hash;
}

// Root cache management functions
static void ClearRootCache() {
    g_rootCache.clear();
}


static void ReserveRootCache(size_t expectedSize) {
    g_rootCache.reserve(expectedSize);
}

// ===== Newton solver implementation =====

// Compute the next starting point based on the polynomial coefficients
static double ComputeStartPoint(const double a[MAX_COEFF], int n) {
    double a0  = log(fabs(a[n]));
    double min = exp((a0 - log(fabs(a[0]))) / static_cast<double>(n));

    for (int i = 1; i < n; i++) {
        if (a[i] != 0.0) {
            double tmp = exp((a0 - log(fabs(a[i]))) / static_cast<double>(n - i));
            if (tmp < min)
                min = tmp;
        }
    }

    double result = min * 0.5;

    return result;
}


// Evaluate a polynomial with real coefficients using Horner's method
// Directly calculates real and imaginary parts without using complex numbers
static Eval EvaluatePolynomial(const double a[MAX_COEFF], int n, double z_re, double z_im) {
    double p = -2.0 * z_re;
    double q = z_re * z_re + z_im * z_im; // Norm of z
    double s = 0.0;
    double r = a[0];
    Eval e;

    for (int i = 1; i < n; i++) {
        double t = a[i] - p * r - q * s;
        s        = r;
        r        = t;
    }

    // Calculate P(z)
    e.z_re  = z_re;
    e.z_im  = z_im;
    e.pz_re = a[n] + z_re * r - q * s;
    e.pz_im = z_im * r;
    e.apz   = sqrt(e.pz_re * e.pz_re + e.pz_im * e.pz_im); // |P(z)|

    return e;
}

// Precomputed constant for machine epsilon (calculate once at startup)
static constexpr double MACHINE_EPSILON = 0.5 * pow((double) _DBL_RADIX, -DBL_MANT_DIG + 1);

// Calculate upper bound for rounding errors (Adam's test)
static double CalculateUpperBound(const double a[MAX_COEFF], int n, double z_re, double z_im) {
    double p = -2.0 * z_re;
    double q = z_re * z_re + z_im * z_im; // Norm of z
    double u = sqrt(q);
    double s = 0.0;
    double r = a[0];
    double e = fabs(r) * (3.5 / 4.5);
    double t;

    for (int i = 1; i < n; i++) {
        t = a[i] - p * r - q * s;
        s = r;
        r = t;
        e = u * e + fabs(t);
    }

    t = a[n] + z_re * r - q * s;
    e = u * e + fabs(t);
    // e = (4.5 * e - 3.5 * (fabs(t) + fabs(r) * u) +
    //      fabs(z_re) * fabs(r)) * 0.5 * pow((double) _DBL_RADIX, -DBL_MANT_DIG + 1);
    e = (4.5 * e - 3.5 * (fabs(t) + fabs(r) * u) +
         fabs(z_re) * fabs(r)) * MACHINE_EPSILON;


    return e;
}

// Real root forward deflation
static void RealDeflation(double a[MAX_COEFF], int &n, double x) {
    double r = 0.0;

    for (int i = 0; i < n; i++)
        a[i]   = r = r * x + a[i];

    n--; // Reduce polynomial degree by 1
}

// Complex root forward deflation for real coefficients
static void ComplexDeflation(double a[MAX_COEFF], int &n, double z_re, double z_im) {
    double r = -2.0 * z_re;
    double u = z_re * z_re + z_im * z_im; // Norm of z

    a[1] -= r * a[0];
    for (int i = 2; i < n - 1; i++)
        a[i]   = a[i] - r * a[i - 1] - u * a[i - 2];

    n -= 2; // Reduce polynomial degree by 2 (complex conjugate roots)
}

// Solve quadratic polynomial - directly write real and imaginary parts
static int SolveQuadratic(const double a[MAX_COEFF], int n,
                          double real_roots[MAX_ROOTS], double imag_roots[MAX_ROOTS],
                          int rootIndex, double new_real[MAX_ROOTS], double new_imag[MAX_ROOTS],
                          int newRootIndex) {
    double real_part, imag_part;
    double r;

    // Notice that a[0] is !=0 since roots=zero has already been handled
    if (n == 1) {
        real_part = -a[1] / a[0];
        imag_part = 0.0;

        real_roots[rootIndex]  = real_part;
        imag_roots[rootIndex]  = imag_part;
        new_real[newRootIndex] = real_part;
        new_imag[newRootIndex] = imag_part;

        return rootIndex + 1;
    } else {
        if (a[1] == 0.0) {
            r = -a[2] / a[0];
            if (r < 0) {
                r         = sqrt(-r);
                real_part = 0.0;
                imag_part = r;

                real_roots[rootIndex]     = real_part;
                imag_roots[rootIndex]     = imag_part;
                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = -imag_part;

                new_real[newRootIndex]     = real_part;
                new_imag[newRootIndex]     = imag_part;
                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = -imag_part;

                return rootIndex + 2;
            } else {
                r = sqrt(r);

                real_roots[rootIndex]     = r;
                imag_roots[rootIndex]     = 0.0;
                real_roots[rootIndex + 1] = -r;
                imag_roots[rootIndex + 1] = 0.0;

                new_real[newRootIndex]     = r;
                new_imag[newRootIndex]     = 0.0;
                new_real[newRootIndex + 1] = -r;
                new_imag[newRootIndex + 1] = 0.0;

                return rootIndex + 2;
            }
        } else {
            r = 1.0 - 4.0 * a[0] * a[2] / (a[1] * a[1]);
            if (r < 0) {
                real_part = -a[1] / (2.0 * a[0]);
                imag_part = a[1] * sqrt(-r) / (2.0 * a[0]);

                real_roots[rootIndex]     = real_part;
                imag_roots[rootIndex]     = imag_part;
                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = -imag_part;

                new_real[newRootIndex]     = real_part;
                new_imag[newRootIndex]     = imag_part;
                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = -imag_part;

                return rootIndex + 2;
            } else {
                real_part = (-1.0 - sqrt(r)) * a[1] / (2.0 * a[0]);
                imag_part = 0.0;

                real_roots[rootIndex] = real_part;
                imag_roots[rootIndex] = imag_part;

                new_real[newRootIndex] = real_part;
                new_imag[newRootIndex] = imag_part;

                real_part = a[2] / (a[0] * real_part);

                real_roots[rootIndex + 1] = real_part;
                imag_roots[rootIndex + 1] = 0.0;

                new_real[newRootIndex + 1] = real_part;
                new_imag[newRootIndex + 1] = 0.0;

                return rootIndex + 2;
            }
        }
    }
}

// C-style interface that directly writes to real and imaginary arrays
// Optimized version without std::complex operations
int PolynomialRootsNewton(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr) {
    int n = degree;      // Polynomial degree
    Eval pz;             // P(z)
    Eval pzprev;         // P(zprev)
    Eval p1z;            // P'(z)
    Eval p1zprev;        // P'(zprev)
    double z_re, z_im;   // Current z (replaces complex<double> z)
    double dz_re, dz_im; // Current step size (replaces complex<double> dz)
    int itercnt;         // Hold the number of iterations per root
    int rootIndex = 0;   // Index for storing roots

    // Constants for complex rotation instead of creating complex numbers
    const double rotation_re = 0.6;
    const double rotation_im = 0.8;

    // Fixed-size array for coefficients
    double coeff[MAX_COEFF];
    memcpy(coeff, coefficient_vector_ptr, sizeof(double) * (degree + 1));

    // Generate hash for this polynomial to check cache - ultra fast integer hash
    uint32_t polyHash     = HashPolynomial(coeff);
    bool usingCachedRoots = false;
    double cachedReal[MAX_ROOTS];
    double cachedImag[MAX_ROOTS];
    int cachedRootCount = 0;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        const RootStorage &storage = cacheIt->second;
        cachedRootCount            = storage.count;
        memcpy(cachedReal, storage.real, sizeof(double) * cachedRootCount);
        memcpy(cachedImag, storage.imag, sizeof(double) * cachedRootCount);
        usingCachedRoots = true;
    }

    // Step 1 eliminate all simple roots
    while (n > 0 && coeff[n] == 0.0) {
        real_zero_vector_ptr[rootIndex]      = 0.0; // Store real part of zero root
        imaginary_zero_vector_ptr[rootIndex] = 0.0; // Store imaginary part of zero root
        rootIndex++;
        n--;
    }

    // Index for tracking cached roots
    int cacheRootIndex = 0;
    double newReal[MAX_ROOTS];
    double newImag[MAX_ROOTS];
    int newRootIndex = 0;

    // Do Newton iteration for polynomial order higher than 2
    while (n > 2) {
        const double Max_stepsize = 5.0; // Allow the next step size to be up to 5x larger
        double r;                        // Current radius
        double rprev;                    // Previous radius
        double eps;                      // Iteration termination value
        bool stage1 = true;              // By default start in stage1
        int steps   = 1;                 // Multisteps if > 1

        // Fixed-size array for derivative coefficients
        double coeffprime[MAX_COEFF - 1];

        // Calculate coefficients of p'(x)
        for (int i        = 0; i < n; i++)
            coeffprime[i] = coeff[i] * double(n - i);

        // Step 2 find a suitable starting point z
        rprev = ComputeStartPoint(coeff, n); // Computed startpoint

        // If we have cached roots, use the next one as a starting point
        if (usingCachedRoots && cacheRootIndex < cachedRootCount) {
            z_re = cachedReal[cacheRootIndex];
            z_im = cachedImag[cacheRootIndex];
            cacheRootIndex++;
        } else {
            // Use default starting point calculation
            if (coeff[n - 1] == 0.0) {
                z_re = 1.0;
                z_im = 0.0;
            } else {
                z_re = -coeff[n] / coeff[n - 1];
                z_im = 0.0;
            }

            // Scale by rprev / |z|
            double abs_z = sqrt(z_re * z_re + z_im * z_im);
            if (abs_z > 0.0) {
                double scale = rprev / abs_z;
                z_re *= scale;
                z_im *= scale;
            }
        }

        // Setup the iteration
        // Current P(z)
        pz = EvaluatePolynomial(coeff, n, z_re, z_im);

        // pzprev which is the previous z or P(0)
        pzprev.z_re  = 0.0;
        pzprev.z_im  = 0.0;
        pzprev.pz_re = coeff[n];
        pzprev.pz_im = 0.0;
        pzprev.apz   = fabs(coeff[n]);

        // p1zprev P'(0) is the P'(0)
        p1zprev.z_re  = pzprev.z_re;
        p1zprev.z_im  = pzprev.z_im;
        p1zprev.pz_re = coeffprime[n - 1];
        p1zprev.pz_im = 0.0;
        p1zprev.apz   = fabs(coeffprime[n - 1]);

        // Set previous dz and calculate the radius of operations.
        dz_re = pz.z_re;
        dz_im = pz.z_im;
        r     = rprev *= Max_stepsize; // Make a reasonable radius of the maximum step size allowed

        // Preliminary eps computed at point P(0) by a crude estimation
        // TODO Epsislon can be pre computed!
        eps = 2 * n * pzprev.apz * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        // Start iteration
        for (itercnt = 0;
             // Check if z + dz != z (equivalent to dz != 0 in complex arithmetic)
             (pz.z_re + dz_re != pz.z_re || pz.z_im + dz_im != pz.z_im) &&
             pz.apz > eps && itercnt < MAX_ITER;
             itercnt++) {
            // Calculate current P'(z)
            p1z = EvaluatePolynomial(coeffprime, n - 1, pz.z_re, pz.z_im);

            // If P'(z)==0 then rotate and try again
            if (p1z.apz == 0.0) {
                // Multiply dz by rotation * Max_stepsize to get away from saddlepoint
                double tmp_re, tmp_im;
                complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                complex_scale(dz_re, dz_im, tmp_re, tmp_im, Max_stepsize);
            } else {
                // Calculate next dz = pz.pz / p1z.pz using complex division
                complex_div(dz_re, dz_im, pz.pz_re, pz.pz_im, p1z.pz_re, p1z.pz_im);

                // Check the Magnitude of Newton's step
                r = sqrt(dz_re * dz_re + dz_im * dz_im);

                if (r > rprev) {
                    // Too large step - rotate and adjust step size
                    double tmp_re, tmp_im, scaled_re, scaled_im;
                    complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                    complex_scale(scaled_re, scaled_im, tmp_re, tmp_im, rprev / r);
                    dz_re = scaled_re;
                    dz_im = scaled_im;
                    r     = sqrt(dz_re * dz_re + dz_im * dz_im);
                }

                // Save 5 times the current step size for the next iteration check
                rprev = r * Max_stepsize;

                // Calculate if stage1 is true or false
                // Stage1 is false if the Newton converges, otherwise true
                double z_re_tmp, z_im_tmp;
                complex_sub(z_re_tmp, z_im_tmp, p1zprev.pz_re, p1zprev.pz_im, p1z.pz_re, p1z.pz_im);
                double denom_re, denom_im;
                complex_sub(denom_re, denom_im, pzprev.z_re, pzprev.z_im, pz.z_re, pz.z_im);

                // z = (p1zprev.pz - p1z.pz) / (pzprev.z - pz.z)
                complex_div(z_re_tmp, z_im_tmp, z_re_tmp, z_im_tmp, denom_re, denom_im);

                // Calculate: (abs(z) / p1z.apz > p1z.apz / pz.apz / 4)
                double abs_z_tmp = sqrt(z_re_tmp * z_re_tmp + z_im_tmp * z_im_tmp);
                stage1           = (abs_z_tmp / p1z.apz > p1z.apz / pz.apz / 4) || (steps != 1);
            }

            // Step accepted. Save pz in pzprev
            pzprev = pz;

            // Next z = pzprev.z - dz
            complex_sub(z_re, z_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);

            // Evaluate polynomial at new z
            pz    = EvaluatePolynomial(coeff, n, z_re, z_im);
            steps = 1;

            if (stage1) {
                // Try multiple steps or shorten steps depending if P(z)
                // is an improvement or not P(z)<P(zprev)
                bool div2;
                double zn_re = pz.z_re;
                double zn_im = pz.z_im;
                Eval npz;

                for (div2 = pz.apz > pzprev.apz; steps <= n; ++steps) {
                    if (div2 == true) {
                        // Shorten steps
                        dz_re *= 0.5;
                        dz_im *= 0.5;
                        complex_sub(zn_re, zn_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);
                    } else {
                        // Try another step in the same direction
                        complex_sub(zn_re, zn_im, zn_re, zn_im, dz_re, dz_im);
                    }

                    // Evaluate new try step
                    npz = EvaluatePolynomial(coeff, n, zn_re, zn_im);
                    if (npz.apz >= pz.apz)
                        break; // Break if no improvement

                    // Improved => accept step and try another round of step
                    pz = npz;

                    if (div2 == true && steps == 2) {
                        // Too many shorten steps => try another direction and break
                        double tmp_re, tmp_im;
                        complex_mul(tmp_re, tmp_im, dz_re, dz_im, rotation_re, rotation_im);
                        dz_re = tmp_re;
                        dz_im = tmp_im;
                        complex_sub(z_re, z_im, pzprev.z_re, pzprev.z_im, dz_re, dz_im);
                        pz = EvaluatePolynomial(coeff, n, z_re, z_im);
                        break;
                    }
                }
            } else {
                // Calculate the upper bound of error using Grant & Hitchins's test
                eps = CalculateUpperBound(coeff, n, pz.z_re, pz.z_im);
            }
        }


        // Check if there is a very small residue in the imaginary part
        // Try to evaluate P(z.real), if that is less than P(z)
        // Take that z.real() is a better root than z
        Eval pztest = EvaluatePolynomial(coeff, n, pz.z_re, 0.0);

        if (pztest.apz <= pz.apz) {
            // Real root
            pz = pztest;

            // Save the root directly to output arrays
            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = 0.0;
            rootIndex++;

            // Also save to the cache arrays
            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = 0.0;
            newRootIndex++;

            RealDeflation(coeff, n, pz.z_re);
        } else {
            // Complex conjugate pair of roots

            // Save both roots directly to output arrays
            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = pz.z_im;
            rootIndex++;

            real_zero_vector_ptr[rootIndex]      = pz.z_re;
            imaginary_zero_vector_ptr[rootIndex] = -pz.z_im;
            rootIndex++;

            // Also save to the cache arrays
            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = pz.z_im;
            newRootIndex++;

            newReal[newRootIndex] = pz.z_re;
            newImag[newRootIndex] = -pz.z_im;
            newRootIndex++;

            ComplexDeflation(coeff, n, pz.z_re, pz.z_im);
        }
    } // End Iteration

    // Solve any remaining linear or quadratic polynomial
    if (n > 0)
        rootIndex = SolveQuadratic(coeff, n, real_zero_vector_ptr, imaginary_zero_vector_ptr,
                                   rootIndex, newReal, newImag, newRootIndex);

    // Store the roots in the cache for future use
    RootStorage storage;
    memcpy(storage.real, newReal, sizeof(double) * newRootIndex);
    memcpy(storage.imag, newImag, sizeof(double) * newRootIndex);
    storage.count         = newRootIndex;
    g_rootCache[polyHash] = storage;

    return rootIndex; // Return the number of roots found
}

// ===== Aberth-Ehrlich solver implementation =====

// Helper function to eliminate zero roots (roots that are exactly zero)
static int ZeroRoots(int n, double a_re[], double a_im[], double res_re[], double res_im[]) {
    int i = n;
    int k = 0;

    // Find roots that are exactly zero (coefficient a[n] = 0)
    while (i > 0 && complex_abs(a_re[i], a_im[i]) < 1e-14) {
        k++;
        i--;
    }

    // If we found zero roots, store them in the result array
    for (int j = n; j > i; j--) {
        res_re[j] = 0.0;
        res_im[j] = 0.0;
    }

    return i; // Return the new degree of the polynomial
}

// Helper function to evaluate a polynomial at a complex point using Horner's method
static Eval Horner(int n, const double a_re[], const double a_im[], double z_re, double z_im) {
    Eval result;
    double result_re = a_re[0];
    double result_im = a_im[0];
    double tmp_re, tmp_im;

    for (int i = 1; i <= n; i++) {
        // result = result * z + a[i]
        complex_mul(tmp_re, tmp_im, result_re, result_im, z_re, z_im);
        complex_add(result_re, result_im, tmp_re, tmp_im, a_re[i], a_im[i]);
    }

    result.z_re  = z_re;
    result.z_im  = z_im;
    result.pz_re = result_re;
    result.pz_im = result_im;
    result.apz   = complex_abs(result_re, result_im);

    return result;
}

// Helper function to calculate starting points for the algorithm
static void StartPoints(int n, const double apolyr[], double Z_re[], double Z_im[]) {
    // Calculate a suitable radius for the starting points
    double r0 = 0.0;
    for (int i = 0; i < n; i++) {
        if (apolyr[i] != 0.0) {
            double tmp = pow(apolyr[i] / apolyr[0], 1.0 / (n - i));
            if (tmp > r0) r0 = tmp;
        }
    }

    // Place starting points uniformly on a circle with radius r0
    const double PI = 3.14159265358979323846;
    r0 *= 1.1; // Slightly larger radius for better convergence

    for (int k = 1; k <= n; k++) {
        double angle = 2 * PI * (k - 1) / n;
        Z_re[k]      = r0 * cos(angle);
        Z_im[k]      = r0 * sin(angle);
    }
}

// Helper function to solve quadratic and linear equations directly
static void Quadratic(int n, double a_re[], double a_im[], double res_re[], double res_im[]) {
    if (n == 1) {
        // Linear equation: a[0]x + a[1] = 0
        complex_div(res_re[1], res_im[1], -a_re[1], -a_im[1], a_re[0], a_im[0]);
    } else if (n == 2) {
        // Quadratic equation: a[0]x^2 + a[1]x + a[2] = 0
        // discriminant = a[1]^2 - 4*a[0]*a[2]
        double disc_re, disc_im;
        double term1_re, term1_im;
        double term2_re, term2_im;
        double sqrt_disc_re, sqrt_disc_im;

        // Calculate a[1]^2
        complex_mul(term1_re, term1_im, a_re[1], a_im[1], a_re[1], a_im[1]);

        // Calculate 4*a[0]*a[2]
        complex_mul(term2_re, term2_im, a_re[0], a_im[0], a_re[2], a_im[2]);
        term2_re *= 4.0;
        term2_im *= 4.0;

        // Calculate discriminant = a[1]^2 - 4*a[0]*a[2]
        complex_sub(disc_re, disc_im, term1_re, term1_im, term2_re, term2_im);

        // Calculate sqrt(discriminant)
        // For complex square root, we use:
        // sqrt(r*e^(i*θ)) = sqrt(r)*e^(i*θ/2)
        double r     = complex_abs(disc_re, disc_im);
        double theta = atan2(disc_im, disc_re);
        sqrt_disc_re = sqrt(r) * cos(theta / 2);
        sqrt_disc_im = sqrt(r) * sin(theta / 2);

        // Calculate -a[1]
        double neg_a1_re = -a_re[1];
        double neg_a1_im = -a_im[1];

        // Calculate 2*a[0]
        double two_a0_re = 2.0 * a_re[0];
        double two_a0_im = 2.0 * a_im[0];

        // Calculate (-a[1] + sqrt(discriminant))/(2*a[0])
        double numer1_re, numer1_im;
        complex_add(numer1_re, numer1_im, neg_a1_re, neg_a1_im, sqrt_disc_re, sqrt_disc_im);
        complex_div(res_re[1], res_im[1], numer1_re, numer1_im, two_a0_re, two_a0_im);

        // Calculate (-a[1] - sqrt(discriminant))/(2*a[0])
        double numer2_re, numer2_im;
        complex_sub(numer2_re, numer2_im, neg_a1_re, neg_a1_im, sqrt_disc_re, sqrt_disc_im);
        complex_div(res_re[2], res_im[2], numer2_re, numer2_im, two_a0_re, two_a0_im);
    }
}

// Main implementation of the Aberth-Ehrlich method
static void AberthEhrlich(int n, const double coeff_re[], const double coeff_im[],
                          double res_re[], double res_im[],
                          const double initialGuesses_re[], const double initialGuesses_im[],
                          int numGuesses) {
    bool dz_flag;
    int itercnt, i, j;
    double f, f0, f1, max_f, eps;

    // Complex variables
    double z_re, z_im, zi_re, zi_im, dz_re, dz_im;

    // Arrays for computation
    double *a_re, *a_im, *w_re, *w_im, *Z_re, *Z_im;
    bool *finish;
    double *apolyr;

    // Temporary variables for complex operations
    double tmp_re, tmp_im;

    // Evaluation results
    Eval fz, fz0, fz1;

    // Oscillation detection variables
    double prev_max_f       = 0.0;
    int oscillation_counter = 0;
    double damping_factor   = 1.0; // Start with no damping

    // Copy the original coefficients
    a_re = new double[n + 1];
    a_im = new double[n + 1];
    for (i = 0; i <= n; i++) {
        a_re[i] = coeff_re[i];
        a_im[i] = coeff_im[i];
    }

    // Eliminate zero roots
    n = ZeroRoots(n, a_re, a_im, res_re, res_im);

    if (n > 2) {
        double *a1_re = new double[n];
        double *a1_im = new double[n];

        // Calculate coefficients of f'(x)
        for (i = 0; i < n; i++) {
            a1_re[i] = a_re[i] * (n - i);
            a1_im[i] = a_im[i] * (n - i);
        }

        w_re   = new double[n + 1];
        w_im   = new double[n + 1];
        apolyr = new double[n + 1];
        Z_re   = new double[n + 1];
        Z_im   = new double[n + 1];
        finish = new bool[n + 1];

        // Simple upper bound for P(z) using Horner with complex coefficients
        f0  = complex_abs(a_re[n], a_im[n]);
        eps = 6 * n * f0 * pow((double) _DBL_RADIX, -DBL_MANT_DIG);

        for (i        = 0; i <= n; i++)
            apolyr[i] = complex_abs(a_re[i], a_im[i]);

        // If we have initial guesses (warm start), use them
        if (initialGuesses_re != nullptr && initialGuesses_im != nullptr && numGuesses >= n) {
            for (i = 1; i <= n; i++) {
                Z_re[i] = initialGuesses_re[i - 1]; // Adjust for 1-indexing in Z
                Z_im[i] = initialGuesses_im[i - 1];
            }
        } else {
            // Otherwise generate starting points
            StartPoints(n, apolyr, Z_re, Z_im);
        }

        for (i        = 1; i <= n; i++)
            finish[i] = false;

        max_f   = 1;
        dz_flag = true;

        // Start iteration
        for (itercnt = 1; dz_flag && max_f > eps && itercnt < MAX_ABERTH_ITER; itercnt++) {
            max_f   = 0;
            dz_flag = false;

            // Early termination for solutions that are "good enough"
            if (max_f < GOOD_ENOUGH_THRESHOLD && itercnt > 3) {
                break;
            }

            // Escalated oscillation handling - if we've been oscillating too long, break out
            if (oscillation_counter > MAX_OSCILLATION_COUNT) {
                break;
            }

            // Check for oscillation
            if (itercnt > 2 && fabs(max_f - prev_max_f) < OSCILLATION_THRESHOLD) {
                oscillation_counter++;
                if (oscillation_counter > 3) {
                    // Apply stronger damping to break oscillation
                    damping_factor *= 0.5;
                    if (damping_factor < 0.01) damping_factor = 0.01; // Minimum damping

                    // Try random perturbation to break symmetry
                    if (oscillation_counter > 4) {
                        for (i = 1; i <= n; i++) {
                            if (!finish[i]) {
                                // Add a tiny random perturbation
                                double scale = 1e-15 * complex_abs(Z_re[i], Z_im[i]);
                                Z_re[i] += scale * ((double) rand() / RAND_MAX - 0.5);
                                Z_im[i] += scale * ((double) rand() / RAND_MAX - 0.5);
                            }
                        }
                    } else {
                    }
                }
            } else {
                // Reset oscillation counter if error is improving significantly
                if (prev_max_f > 0 && max_f < 0.5 * prev_max_f) {
                    oscillation_counter = 0;
                    // Gradually restore damping if not oscillating
                    if (damping_factor < 1.0) damping_factor *= 1.2;
                    if (damping_factor > 1.0) damping_factor = 1.0;
                }
            }

            prev_max_f = max_f;

            // Count how many roots are left to find
            int unfinished_roots = 0;
            for (i = 1; i <= n; i++) {
                if (!finish[i]) unfinished_roots++;
            }

            // If we're down to just a few roots, check if we can finish them directly
            if (unfinished_roots <= 2 && itercnt > 10) {
                // Try to solve directly for the remaining roots
                bool all_done = true;
                for (i = 1; i <= n; i++) {
                    if (!finish[i]) {
                        zi_re = Z_re[i];
                        zi_im = Z_im[i];
                        fz0   = Horner(n, a_re, a_im, zi_re, zi_im);
                        if (fz0.apz < 1e-10) {
                            finish[i] = true;
                        } else {
                            all_done = false;
                        }
                    }
                }
                if (all_done) break;
            }

            for (i = 1; i <= n; i++) {
                if (finish[i] == true) continue;

                zi_re = Z_re[i];
                zi_im = Z_im[i];
                fz0   = Horner(n, a_re, a_im, zi_re, zi_im);
                f0    = fz0.apz;
                fz1   = Horner(n - 1, a1_re, a1_im, zi_re, zi_im);
                f1    = fz1.apz;

                // Calculate w[i] = sum(1 / (zi - Z[j])) for j != i
                w_re[i] = 0.0;
                w_im[i] = 0.0;

                for (j = 1; j <= n; j++) {
                    if (i != j) {
                        // dz = zi - Z[j]
                        complex_sub(dz_re, dz_im, zi_re, zi_im, Z_re[j], Z_im[j]);

                        // tmp = 1 / dz
                        complex_div(tmp_re, tmp_im, 1.0, 0.0, dz_re, dz_im);

                        // w[i] += tmp
                        complex_add(w_re[i], w_im[i], w_re[i], w_im[i], tmp_re, tmp_im);
                    }
                }

                // dz = fz1 / fz0 - w[i]
                complex_div(tmp_re, tmp_im, fz1.pz_re, fz1.pz_im, fz0.pz_re, fz0.pz_im);
                complex_sub(dz_re, dz_im, tmp_re, tmp_im, w_re[i], w_im[i]);

                // dz = 1 / dz
                complex_div(dz_re, dz_im, 1.0, 0.0, dz_re, dz_im);

                // Apply damping factor to prevent oscillation
                dz_re *= damping_factor;
                dz_im *= damping_factor;

                // Store dz in w[i] for later use
                w_re[i] = dz_re;
                w_im[i] = dz_im;

                // z = zi - dz
                complex_sub(z_re, z_im, zi_re, zi_im, dz_re, dz_im);

                // Evaluate at new point
                fz = Horner(n, a_re, a_im, z_re, z_im);
                f  = fz.apz;

                // Add line search to ensure improvement
                if (f > f0) {
                    // If the step increases error, try more step sizes
                    double best_alpha = 1.0;
                    double best_f     = f;
                    double best_z_re  = z_re;
                    double best_z_im  = z_im;

                    // Try more step sizes (0.8, 0.6, 0.4, 0.2, 0.1, 0.01)
                    const double alphas[] = {0.8, 0.6, 0.4, 0.2, 0.1, 0.01};
                    for (const double alpha: alphas) {
                        // z_new = zi - alpha * dz
                        double z_new_re, z_new_im;
                        complex_scale(tmp_re, tmp_im, dz_re, dz_im, alpha);
                        complex_sub(z_new_re, z_new_im, zi_re, zi_im, tmp_re, tmp_im);

                        // Evaluate at new point
                        Eval fz_new  = Horner(n, a_re, a_im, z_new_re, z_new_im);
                        double f_new = fz_new.apz;

                        if (f_new < best_f) {
                            best_f     = f_new;
                            best_z_re  = z_new_re;
                            best_z_im  = z_new_im;
                            best_alpha = alpha;
                        }
                    }

                    // Use the best step size
                    z_re = best_z_re;
                    z_im = best_z_im;
                    f    = best_f;

                    // If we're taking a very small step, that's another sign of oscillation
                    if (best_alpha <= 0.1) {
                        oscillation_counter++;
                    }
                }

                Z_re[i] = z_re;
                Z_im[i] = z_im;

                // Check if dz is effectively zero (z + dz == z)
                dz_flag = dz_flag || (z_re + dz_re != z_re || z_im + dz_im != z_im);

                if (f > max_f)
                    max_f = f;

                if (f <= eps || (z_re + dz_re == z_re && z_im + dz_im == z_im)) {
                    double z0_re, z0_im;
                    finish[i] = true;

                    if (fabs(z_re) >= fabs(z_im)) {
                        z0_re = z_re;
                        z0_im = 0.0;
                    } else {
                        z0_re = 0.0;
                        z0_im = z_im;
                    }

                    Eval fz0 = Horner(n, a_re, a_im, z0_re, z0_im);

                    if (fz0.apz <= f) {
                        Z_re[i] = z_re = z0_re;
                        Z_im[i] = z_im = z0_im;
                    }
                }
            }
        }

        for (i = 1; i <= n; i++) {
            res_re[i] = Z_re[i];
            res_im[i] = Z_im[i];
        }

        delete[] finish;
        delete[] Z_re;
        delete[] Z_im;
        delete[] w_re;
        delete[] w_im;
        delete[] a1_re;
        delete[] a1_im;
        delete[] apolyr;
    } else {
        Quadratic(n, a_re, a_im, res_re, res_im);
    }

    delete[] a_re;
    delete[] a_im;
    return;
}

// C-style interface that directly writes to real and imaginary arrays
// Designed to match the interface of the Newton solver
int PolynomialRootsAberth(
    const double *coefficient_vector_ptr,
    int degree,
    double *real_zero_vector_ptr,
    double *imaginary_zero_vector_ptr) {
    // Generate hash for this polynomial to check cache
    uint32_t polyHash     = HashPolynomial(coefficient_vector_ptr);
    bool usingCachedRoots = false;
    double cachedReal[MAX_ROOTS];
    double cachedImag[MAX_ROOTS];
    int cachedRootCount = 0;

    // Check if we have solved a similar polynomial before
    auto cacheIt = g_rootCache.find(polyHash);
    if (cacheIt != g_rootCache.end()) {
        const RootStorage &storage = cacheIt->second;
        cachedRootCount            = storage.count;
        memcpy(cachedReal, storage.real, sizeof(double) * cachedRootCount);
        memcpy(cachedImag, storage.imag, sizeof(double) * cachedRootCount);
        usingCachedRoots = true;
    }

    // Convert real coefficients to separate real/imag arrays
    double complexCoeffs_re[MAX_COEFF];
    double complexCoeffs_im[MAX_COEFF];
    for (int i = 0; i <= degree; i++) {
        complexCoeffs_re[i] = coefficient_vector_ptr[i];
        complexCoeffs_im[i] = 0.0;
    }

    // Aberth implementation is 1-indexed
    double roots_re[MAX_ROOTS + 1];
    double roots_im[MAX_ROOTS + 1];

    // Call the Aberth-Ehrlich implementation, with warm start if available
    if (usingCachedRoots) {
        AberthEhrlich(degree, complexCoeffs_re, complexCoeffs_im, roots_re, roots_im,
                      cachedReal, cachedImag, cachedRootCount);
    } else {
        AberthEhrlich(degree, complexCoeffs_re, complexCoeffs_im, roots_re, roots_im);
    }

    // Convert results to the output arrays
    int rootCount = 0;
    double newReal[MAX_ROOTS];
    double newImag[MAX_ROOTS];

    for (int i = 1; i <= degree; i++) {
        real_zero_vector_ptr[rootCount]      = roots_re[i];
        imaginary_zero_vector_ptr[rootCount] = roots_im[i];

        newReal[rootCount] = roots_re[i];
        newImag[rootCount] = roots_im[i];

        rootCount++;
    }

    // Store the roots in the cache for future use
    RootStorage storage;
    memcpy(storage.real, newReal, sizeof(double) * rootCount);
    memcpy(storage.imag, newImag, sizeof(double) * rootCount);
    storage.count         = rootCount;
    g_rootCache[polyHash] = storage;

    return rootCount;
}
