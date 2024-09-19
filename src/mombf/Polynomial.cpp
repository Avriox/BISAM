//=======================================================================
// Copyright (C) 2003-2013 William Hallahan
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without restriction,
// including without limitation the rights to use, copy, modify, merge,
// publish, distribute, sublicense, and/or sell copies of the Software,
// and to permit persons to whom the Software is furnished to do so,
// subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
// OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
// HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
// WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
// FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
// OTHER DEALINGS IN THE SOFTWARE.
//=======================================================================

//**********************************************************************
//  File: Polynomial.cpp
//  Author: Bill Hallahan
//  Date: January 30, 2003
//
//  Abstract:
//
//    This file contains the implementation for class Polynomial.
//
//**********************************************************************


#include "Polynomial.h"
#include "PolynomialRootFinder.h"

//======================================================================
//  Constructor: Polynomial::Polynomial
//======================================================================

Polynomial::Polynomial()
        : m_degree(-1), m_coefficient_vector_ptr(NULL) {
    SetToScalar(0.0);
}


//======================================================================
//  Destructor: Polynomial::~Polynomial
//======================================================================

Polynomial::~Polynomial() {
}

//======================================================================
//  Member Function: Polynomial::SetCoefficients
//
//  Abstract:
//
//    This method sets the polynomial coefficients.
//
//
//  Input:
//
//    coefficient_vector_ptr    The vector of coefficients in order
//                              of increasing powers.
//
//    degree                    The degree of the polynomial.
//
//
//  Return Value:
//
//    The function has no return value.
//
//======================================================================

void Polynomial::SetCoefficients(double *coefficient_vector_ptr,
                                 int degree) {

    m_degree = degree;

    SetLength(m_degree + 1, false);

    int ii = 0;

    for (ii = 0; ii <= m_degree; ++ii) {
        m_coefficient_vector_ptr[ii] = coefficient_vector_ptr[ii];
    }

    AdjustPolynomialDegree();
}

//======================================================================
//  Member Function: Polynomial::SetToScalar
//
//  Abstract:
//
//    This method sets the polynomial to be a scalar.
//    The polynomial degree is set to 0 in this method.
//
//
//  Input:
//
//    scalar                    A scalar value
//
//  Return Value:
//
//    The function has no return value.
//
//======================================================================

void Polynomial::SetToScalar(double scalar) {
    SetCoefficients(&scalar, 0);
}


//======================================================================
//  Member Function: Polynomial::FindRoots
//
//  Abstract:
//
//    This method determines the roots of a polynomial which has
//    real coefficients.
//
//
//  Input:
//
//
//    real_zero_vector_ptr       A double precision vector that will
//                               contain the real parts of the roots
//                               of the polynomial when this method
//                               returns.
//
//    imaginary_zero_vector_ptr  A double precision vector that will
//                               contain the real parts of the roots
//                               of the polynomial when this method
//                               returns.
//
//    roots_found_ptr           A pointer to an integer that will
//                              equal the number of roots found when
//                              this method returns. If the method
//                              returns SUCCESS then this value will
//                              always equal the degree of the
//                              polynomial.
//
//  Return Value:
//
//    This function returns an enum value of type
//    'PolynomialRootFinder::RootStatus_T'.
//
//======================================================================

PolynomialRootFinder::RootStatus_T Polynomial::FindRoots(double *real_zero_vector_ptr,
                                                         double *imaginary_zero_vector_ptr,
                                                         int *roots_found_ptr) const {
    //assert(m_degree >= 0);

    PolynomialRootFinder *polynomial_root_finder_ptr = new PolynomialRootFinder;

    if (polynomial_root_finder_ptr == NULL) {
        throw std::bad_alloc();
    }

    std::unique_ptr<PolynomialRootFinder> root_finder_ptr(polynomial_root_finder_ptr);
    //std::auto_ptr<PolynomialRootFinder> root_finder_ptr(polynomial_root_finder_ptr);  //auto_ptr deprecated and replaced by unique_ptr

    PolynomialRootFinder::RootStatus_T status = root_finder_ptr->FindRoots(m_coefficient_vector_ptr,
                                                                           m_degree,
                                                                           real_zero_vector_ptr,
                                                                           imaginary_zero_vector_ptr,
                                                                           roots_found_ptr);
    return status;
}


//======================================================================
//  Member Function: Polynomial::AdjustPolynomialDegree
//
//  Abstract:
//
//    This method will decrease the polynomial degree until leading
//    coefficient is non-zero or until the polynomial degree is zero.
//
//
//  Input:
//
//    None.
//
//
//  Return Value:
//
//    This method has no return value.
//
//======================================================================

void Polynomial::AdjustPolynomialDegree() {
    //------------------------------------------------------------------
    //  Any leading coefficient with a magnitude less than DBL_EPSILON
    //  is treated as if it was zero.
    //------------------------------------------------------------------

    while ((m_degree > 0)
           && (fabs(m_coefficient_vector_ptr[m_degree]) < DBL_EPSILON)) {
        m_coefficient_vector_ptr[m_degree] = 0.0;
        m_degree--;
    }

    return;
}


//======================================================================
//  Member Function: Polynomial::SetLength
//
//  Abstract:
//
//    This function is called to set the buffer lengths for the
//    coefficient vector. If the new buffer length is less than
//    or equal to the current buffer lengths then then the buffer
//    is not modified and the data length is set to the new buffer
//    length. If the new data length is greater than the current
//    buffer lengths then the buffer is reallocated to the new
//    buffer size. In this case, if argument copy_data_flag
//    is set to the value true then the data in the old buffer
//    is copied to the new buffer.
//
//
//  Input:
//
//    udata_length             The new length of the data.
//
//    copy_data_flag           If this is true then existing data
//                             is copied to any newly allocated buffer.
//                             This parameter defaults to the value
//                             'true'.
//
//  Output:
//
//    This function has no return value.
//
//======================================================================

void Polynomial::SetLength(unsigned int number_of_coefficients,
                           bool copy_data_flag) {

    // If m_degree is equal to -1, then this is a new polynomial and the
    // caller will set m_degree.
    if ((!copy_data_flag) || (m_degree == -1)) {
        // Clear and resize the coefficient vector.
        m_coefficient_vector.clear();
        m_coefficient_vector.resize(number_of_coefficients);
        m_coefficient_vector_ptr = &m_coefficient_vector[0];
    } else {
        // Save the polynomial values in a temporary vector.
        std::vector<double> temp_vector;
        temp_vector.resize(m_degree + 1);

        int i = 0;

        for (i = 0; i <= m_degree; ++i) {
            temp_vector[i] = m_coefficient_vector_ptr[i];
        }

        // Clear and resize the coefficient vector.
        m_coefficient_vector.clear();
        m_coefficient_vector.resize(number_of_coefficients);
        m_coefficient_vector_ptr = &m_coefficient_vector[0];

        // Restore the coefficients for the new vector size.
        // Was the polynomial size increased?
        if (number_of_coefficients > (unsigned int) (m_degree + 1)) {
            // The polynomial size was increased.
            for (i = 0; i <= m_degree; ++i) {
                m_coefficient_vector_ptr[i] = temp_vector[i];
            }

            for (i = m_degree + 1; i < (int) (number_of_coefficients); ++i) {
                m_coefficient_vector_ptr[i] = 0.0;
            }
        } else {
            for (int i = 0; i < (int) (number_of_coefficients); ++i) {
                m_coefficient_vector_ptr[i] = temp_vector[i];
            }
        }
    }

    return;
}

