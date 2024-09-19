/***********************************************************
 Basic, input-output and matrix manipulation

 Authors: Peter Mueller, Stephen Morris, David Rossell
          (some routines obtained from other sources)
 Edits: P. Roebuck
***********************************************************/

#include <iostream>
#include "cstat.h"



/*
 * Much of this code has undocumented assumptions, the least of which
 * is IEC-559 / IEEE-754 standard compliance.
 */


/*
 * Globals
 */
static long is1 = 123456789, is2 = 981963;
static int cstat_set = 0;

FILE *ifile, *ofile;
int nv = 0;

long Xm1, Xm2;
long Xa1, Xa2;
long Xcg1[32], Xcg2[32];
long Xa1w, Xa2w;
long Xig1[32], Xig2[32];
long Xlg1[32], Xlg2[32];
long Xa1vw, Xa2vw;
long Xqanti[32];



/******************************************************************************
                          MEMORY ALLOCATION
******************************************************************************/

/* Allocate int vector with subscript range v[nl..nh] */
int *ivector(int nl,
             int nh) {
    int *v;
    size_t count = nh - nl + 1;

    //assert(count >= 0);

    nv += count;
    v = (int *) calloc(count, sizeof(int));
    if (v == NULL) {
//        nrerror("ivector", "allocate an int vector", "");
        /*NOTREACHED*/
    }
    return v - nl;
}


/* Allocate double vector with subscript range v[nl..nh] */
double *dvector(int nl,
                int nh) {
    double *v;
    size_t count = nh - nl + 1;

    //assert(count >= 0);

    nv += count;
    v = (double *) calloc(count, sizeof(double));
    if (v == NULL) {
//        nrerror("dvector", "allocate a double vector", "");
        /*NOTREACHED*/
    }
    return v - nl;
}


/* Allocate double matrix with subscript range m[nrl..nrh][ncl..nch] */
double **dmatrix(int nrl, int nrh, int ncl, int nch) {
    double **m;
    size_t nrow = nrh - nrl + 1;
    size_t ncol = nch - ncl + 1;
    int i;

    nv += nrow * ncol;

    /* Allocate pointers to rows */
    m = (double **) calloc(nrow, sizeof(double *));
    if (m == NULL) {
//        nrerror("dmatrix", "allocate a double matrix (1st dim)", "");
        /*NOTREACHED*/
    }
    m -= nrl;

    /* For each row pointer... */
    for (i = nrl; i <= nrh; i++) {
        /* Allocate columns for individual row */
        m[i] = (double *) calloc(ncol, sizeof(double));
        if (m[i] == NULL) {
//            nrerror("dmatrix", "allocate a double matrix (2nd dim)", "");
            /*NOTREACHED*/
        }
        m[i] -= ncl;
    }
    return m;
}


/* Free int vector allocated with ivector() */
void free_ivector(int *v,
                  int nl,
                  int nh) {

    //if ((v+nl) != NULL) {
    free((char *) (v + nl));
    //}
    nv -= (nh - nl + 1);
}


/* Free double vector allocated with dvector() */
void free_dvector(double *v, int nl, int nh) {

    //if ((v+nl) != NULL) {
    free((char *) (v + nl));
    //}
    nv -= (nh - nl + 1);
}


/* Free double matrix allocated by dmatrix() */
void free_dmatrix(double **m, int nrl, int nrh, int ncl, int nch) {
    int i;
    size_t nrow = nrh - nrl + 1;
    size_t ncol = nch - ncl + 1;

    for (i = nrh; i >= nrl; i--) {
        //if ((m[i]+ncl) != NULL) {
        free((char *) (m[i] + ncl));
        //}
    }
    //if ((m+nrl) != NULL) {
    free((char *) (m + nrl));
    //}
    nv -= ncol * nrow;
}



/************************************************************************
                          MATHEMATICAL FUNCTIONS
************************************************************************/



/*
-----------------------------------------------------------------------
            EVALUATION OF LN(GAMMA(A)) FOR POSITIVE A
-----------------------------------------------------------------------
     WRITTEN BY ALFRED H. MORRIS
          NAVAL SURFACE WARFARE CENTER
          DAHLGREN, VIRGINIA
--------------------------
     D = 0.5*(LN(2*PI) - 1)
--------------------------
*/
double gamln(double *a) {
    static double c0 = .833333333333333e-01;
    static double c1 = -.277777777760991e-02;
    static double c2 = .793650666825390e-03;
    static double c3 = -.595202931351870e-03;
    static double c4 = .837308034031215e-03;
    static double c5 = -.165322962780713e-02;
    static double d = .418938533204673e0;
    static double gamln, t, w;
    static int i, n;
    static double T1;
/*
     ..
     .. Executable Statements ..
*/
    if (*a > 0.8e0) goto S10;
    gamln = gamln1(a) - log(*a);
    return gamln;
    S10:
    if (*a > 2.25e0) goto S20;
    t = *a - 0.5e0 - 0.5e0;
    gamln = gamln1(&t);
    return gamln;
    S20:
    if (*a >= 10.0e0) goto S40;
    n = (int) (*a - 1.25e0);
    //n = (long)(*a - 1.25e0);
    t = *a;
    w = 1.0e0;
    for (i = 1; i <= n; i++) {
        t -= 1.0e0;
        w = t * w;
    }
    T1 = t - 1.0e0;
    gamln = gamln1(&T1) + log(w);
    return gamln;
    S40:
    t = pow(1.0e0 / *a, 2.0);
    w = (((((c5 * t + c4) * t + c3) * t + c2) * t + c1) * t + c0) / *a;
    gamln = d + w + (*a - 0.5e0) * (log(*a) - 1.0e0);
    return gamln;
}


/*
-----------------------------------------------------------------------
     EVALUATION OF LN(GAMMA(1 + A)) FOR -0.2 .LE. A .LE. 1.25
-----------------------------------------------------------------------
*/
double gamln1(double *a) {
    static double p0 = .577215664901533e+00;
    static double p1 = .844203922187225e+00;
    static double p2 = -.168860593646662e+00;
    static double p3 = -.780427615533591e+00;
    static double p4 = -.402055799310489e+00;
    static double p5 = -.673562214325671e-01;
    static double p6 = -.271935708322958e-02;
    static double q1 = .288743195473681e+01;
    static double q2 = .312755088914843e+01;
    static double q3 = .156875193295039e+01;
    static double q4 = .361951990101499e+00;
    static double q5 = .325038868253937e-01;
    static double q6 = .667465618796164e-03;
    static double r0 = .422784335098467e+00;
    static double r1 = .848044614534529e+00;
    static double r2 = .565221050691933e+00;
    static double r3 = .156513060486551e+00;
    static double r4 = .170502484022650e-01;
    static double r5 = .497958207639485e-03;
    static double s1 = .124313399877507e+01;
    static double s2 = .548042109832463e+00;
    static double s3 = .101552187439830e+00;
    static double s4 = .713309612391000e-02;
    static double s5 = .116165475989616e-03;
    static double gamln1, w, x;
/*
     ..
     .. Executable Statements ..
*/
    if (*a >= 0.6e0) goto S10;
    w = ((((((p6 * *a + p5) * *a + p4) * *a + p3) * *a + p2) * *a + p1) * *a + p0) / ((((((q6 * *a + q5) * *a +
                                                                                          q4) * *a + q3) * *a + q2) *
                                                                                       *a + q1) * *a + 1.0e0);
    gamln1 = -(*a * w);
    return gamln1;
    S10:
    x = *a - 0.5e0 - 0.5e0;
    w = (((((r5 * x + r4) * x + r3) * x + r2) * x + r1) * x + r0) / (((((s5 * x + s4) * x + s3) * x + s2) * x + s1) * x
                                                                     + 1.0e0);
    gamln1 = x * w;
    return gamln1;
}


/* Bernoulli numbers of even order from 2 to 60 */
static double bernou[30] = {
        1.0 / 6.0,
        -1.0 / 30.0,
        1.0 / 42.0,
        -1.0 / 30.0,
        5.0 / 66.0,
        -691.0 / 2730.0,
        7.0 / 6.0,
        -3617.0 / 510.0,
        43867.0 / 798.0,
        -174611.0 / 330.0,
        854513.0 / 138.0,
        -236364091.0 / 2730.0,
        8553103.0 / 6.0,
        -23749461029.0 / 870.0,
        8615841276005.0 / 14322.0,
        -7709321041217.0 / 510.0,
        2577687858367.0 / 6.0,
        -1.371165521e13,
        4.883323190e14,
        -1.929657934e16,
        8.416930476e17,
        -4.033807185e19,
        2.115074864e21,
        -1.208662652e23,
        7.500866746e24,
        -5.038778101e26,
        3.652877648e28,
        -2.849876930e30,
        2.386542750e32,
        -2.139994926e34
};


/* log of Beta function */
double lnbeta(double a,
              double b) {
    double c = a + b;
    return (gamln(&a) + gamln(&b) - gamln(&c));
}


///* Returns 1.0 if x>=0, -1.0 if x<0 */
double dsign(double x) {
    return (x >= 0) ? 1.0 : -1.0;
}

///************************************************************************
//                            VECTOR ALGEBRA
//************************************************************************/


/*
 * Multiply matrix A[ini..fi][ini..fi] by vector x[ini..fi]
 * and add vector y[ini..fi].
 * Store result in vector z.
 */
void Ax_plus_y(double **A,
               const double *x,
               const double *y,
               double *z,
               int ini,
               int fi) {
    int i, j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(y != NULL);
    //assert(z != NULL);

    for (i = ini; i <= fi; i++) {
        z[i] = y[i];
        for (j = ini; j <= fi; j++) {
            z[i] += A[i][j] * x[j];
        }
    }
}


/*
 * Multiply matrix A[rowini..rowfi][colini..colfi] by vector x[colini..colfi].
 * Store result in vector z.
 */
void Ax(double **A,
        const double *x,
        double *z,
        int rowini,
        int rowfi,
        int colini,
        int colfi) {
    int i;
    int j;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[i][j] * x[j];
        }
    }
}


void Avecx_int(const int *A, const double *x, double *z, int rowini, int rowfi, int colini, int colfi) {
    int i;
    int j;
    int nrow = rowfi - rowini + 1;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[i + j * nrow] * x[j];
        }
    }
}


void Atvecx_int(const int *A,
                const double *x,
                double *z,
                int rowini,
                int rowfi,
                int colini,
                int colfi) {
    int i;
    int j;
    int ncol = colfi - colini + 1;

    //assert(A != NULL);
    //assert(x != NULL);
    //assert(z != NULL);

    for (i = rowini; i <= rowfi; i++) {
        z[i] = 0.0;
        for (j = colini; j <= colfi; j++) {
            z[i] += A[j + i * ncol] * x[j];
        }
    }
}


/*
 * Returns sum of multiplying symmetric (implicit matrix) vector A[sel][sel]
 * by transposed vector x[0..nsel] by vector x[0..nsel] for quadratic forms.
 *     ncolA: number of columns in A.
 *     nsel : length of vector sel.
 *     sel  : vector with indexes for (rows,cols) in A to be used in the operation.
 *
 * Same as above but subset is only for A.
 * Note: Faster than xtAy() for symmetric A (saves 25%-50% operations).
 */
double quadratic_xtAselx(const double *x,
                         crossprodmat *A,
                         const int *ncolA,
                         const int *nsel,
                         const int *sel) {
    int i;
    int j;
    double z = 0.0;

    //assert(x != NULL);
    //assert(A != NULL);
    //assert(ncolA != NULL);
    //assert(nsel != NULL);
    //assert(sel != NULL);

    for (i = 0; i <= (*nsel) - 1; i++) {
        int i_sel;

        i_sel = sel[i];
        z += (A->at(i_sel * (*ncolA) + i_sel)) * x[i] * x[i];
        for (j = i + 1; j <= (*nsel) - 1; j++) {
            z += 2 * (A->at(i_sel * (*ncolA) + sel[j])) * x[i] * x[j];
        }
    }
    return (z);
}

// Implementation of quadratic_xtAselx using Eigen
double quadratic_xtAselx_eigen(const Eigen::VectorXd &x, const Eigen::MatrixXd &A, const Eigen::VectorXi &sel) {
    Eigen::VectorXd x_sel = x;
    Eigen::MatrixXd A_sel = A(sel, sel);
    return x_sel.transpose() * A_sel * x_sel;
}


void AvectBvec_int(int *A, int nrowA, int ncolA, int *B, int nrowB, int ncolB, double **C) {
    int i;
    int j;
    int k;

    //assert(A != NULL);
    //assert(B != NULL);
    //assert(C != NULL);

    if (nrowA != nrowB) {
//        errorC("AvectBvec", "dimensions don't match", 1);
        /*NOTREACHED*/
    }
    for (i = 1; i <= ncolA; i++) {
        int offsetA = (i - 1) * nrowA;
        for (j = 1; j <= ncolB; j++) {
            C[i][j] = 0.0;
            int offsetB = (j - 1) * nrowB;
            for (k = 1; k <= nrowA; k++) {
                C[i][j] += A[k - 1 + offsetA] * B[k - 1 + offsetB];
            }
        }
    }
}


///* Diagonal matrix */
void ddiag(double **A,
           int ini,
           int fi) {
    int i;
    int j;

    //assert(A != NULL);

    for (i = ini; i <= fi; i++) {
        for (j = ini; j <= fi; j++) {
            A[i][j] = (i == j) ? 1 : 0;
        }
    }
}


// TODO WHY THE FUCK WOULD YOU WRITE YOUR OWN VERISON
double max_xy(double x,
              double y) {
    return (x > y) ? x : y;
}


void choldc(double **a, int n, double **aout, bool *posdef) {
/* Given a positive-definite symmetric matrix a[1..n][1..n], this routine constructs its Cholesky
decomposition, A = L * L' . On input, only the upper triangle of a need be given;
 The Cholesky factor L is returned in the lower triangle of aout (upper-diag elem are set to 0) */
    int i, j, k;
    double sum, *p, max_a;

    *posdef = true;
    for (i = 1; i <= n; i++) { for (j = i; j <= n; j++) { aout[i][j] = a[i][j]; }}  //copy a into aout
    p = dvector(1, n);
    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            for (sum = aout[i][j], k = i - 1; k >= 1; k--) sum -= aout[i][k] * aout[j][k];
            if (i == j) {
                if (sum <= 0.0) *posdef = false;
                aout[i][i] = sqrt(sum);
            } else {
                max_a = max_xy(fabs(aout[i][i]), 1e-10);
                aout[j][i] = sum / max_a;
            }
        }
    }
    free_dvector(p, 1, n);
    for (i = 1; i <= n; i++) { for (j = i + 1; j <= n; j++) { aout[i][j] = 0; }}  //set upper-diagonal elem to 0
}


void choldc_inv(double **a, int n, double **aout, bool *posdef) {
    /*Given a positive-definite symmetric matrix a[1..n][1..n], this routine computes the inverse
   of its Cholesky matrix. That is, if A=L * L' it returns the inverse of L
   (note that inv(A)= inv(L)' * inv(L)) */
    choldc(a, n, aout, posdef);
    if (*posdef) {
        choldc_inv_internal(aout, n);
    }
}


void cholS_inv(double **cholS, int n, double **cholSinv) {
    /*Given the Cholesky decomposition of a matrix S, which we denote cholS, returns the inverse of cholS */
    int i, j;
    for (i = 1; i <= n; i++) for (j = 1; j <= i; j++) cholSinv[i][j] = cholS[i][j];
    choldc_inv_internal(cholSinv, n);
}

void choldc_inv_internal(double **cholS, int n) {
    /*Computes inverse of Cholesky matrix cholS and stores the result in cholS*/
    int i, j, k;
    double sum, max_a;
    for (i = 1; i <= n; i++) {
        max_a = max_xy(cholS[i][i], 1e-10);
        cholS[i][i] = 1.0 / max_a;
        for (j = i + 1; j <= n; j++) {
            sum = 0.0;
            for (k = i; k < j; k++) sum -= cholS[j][k] * cholS[k][i];
            max_a = max_xy(cholS[j][j], 1e-10);
            cholS[j][i] = sum / max_a;
        }
    }
}


/*
 * Find determinant of the matrix having chols as its Cholesky decomposition.
 *
 * Example:
 *   choldc(S, n, cholS, posdef);
 *   det = choldc_det(cholS, n);
 *
 * Another example:
 *   choldc_inv(S, n, cholSinv, posdef);
 *   det = 1.0 / choldc_det(cholSinv, n);
 */
double choldc_det(double **chols, int n) {
    int i;
    double value, det = 1.0;

    //assert(chols != NULL);
    for (i = 1; i <= n; i++) {
        value = chols[i][i];
        det *= value * value;
    }
    return (det);
}

double logcholdc_det(double **chols, int n) {
    int i;
    double logdet = 0;

    //assert(chols != NULL);
    for (i = 1; i <= n; i++) { logdet += log(chols[i][i]); }
    return (2.0 * logdet);
}


/*
  Inverse of a symmetric, positive definite matrix a[1..n][1..n] using Cholesky decomposition.

  Input: either a, its Cholesky decomposition chola, or the inverse of its Cholesky decomposition cholainv. If chola is provided then a is ignored. If cholainv is provided, then chola is ignored.

  Output: aout contains the inverse of a, posdef returns if matrix was indeed positive definite

 */
void inv_posdef(double **a, int n, double **aout, bool *posdef, double **chola, double **cholainv) {
    int i, j;
    double **b;

    if (cholainv == NULL) {
        b = dmatrix(1, n, 1, n);
        if (chola == NULL) {
            choldc_inv(a, n, b, posdef); //inverse of chol(a)
        } else {
            cholS_inv(chola, n, b);  //inverse of chola
        }
    } else {
        b = cholainv;
    }

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) { sum += b[k][i] * b[k][j]; }
            aout[i][j] = sum;
        }
    }

    if (cholainv == NULL) free_dmatrix(b, 1, n, 1, n);

    for (i = 2; i <= n; i++) {
        for (j = 1; j < i; j++) {
            aout[i][j] = aout[j][i];
        }
    }
}


/*
 * Inverse of a symmetric, positive definite matrix a[1..n][1..n] using
 * Cholesky decomposition. Result is returned in aout.
 * Does the same as inv_posdef, except that here only upper triangular
 * elements are returned.
 */
void inv_posdef_upper(double **a,
                      int n,
                      double **aout,
                      bool *posdef) {
    int i;
    int j;
    double **b;

    //assert(a != NULL);
    //assert(aout != NULL);

    b = dmatrix(1, n, 1, n);
    choldc_inv(a, n, b, posdef);
    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) {
                sum += b[k][i] * b[k][j];
            }
            aout[i][j] = sum;
        }
    }
    free_dmatrix(b, 1, n, 1, n);
}


/*
 * Inverse and determinant of a positive definite matrix a[1..n][1..n] using
 * Cholesky decomposition. Inverse is returned in aout, determinant in det_a.
 */
void invdet_posdef(double **a,
                   int n,
                   double **aout,
                   double *det_a) {
    bool posdef;
    int i;
    int j;
    double **b;

    //assert(a != NULL);
    //assert(aout != NULL);
    //assert(det_a != NULL);

    b = dmatrix(1, n, 1, n);
    choldc_inv(a, n, b, &posdef);
    *det_a = 1.0;
    for (i = 1; i <= n; i++) {
        double value;

        value = b[i][i];
        (*det_a) *= 1 / (value * value);
    }

    for (i = 1; i <= n; i++) {
        for (j = i; j <= n; j++) {
            int k;
            double sum;

            sum = 0.0;
            for (k = 1; k <= n; k++) {
                sum += b[k][i] * b[k][j];
            }
            aout[i][j] = sum;
        }
    }
    free_dmatrix(b, 1, n, 1, n);

    for (i = 2; i <= n; i++) {
        for (j = 1; j < i; j++) {
            aout[i][j] = aout[j][i];
        }
    }
}


/*
 * Sorts vector of doubles x by rearranging values in index with quicksort
 * algorithm, e.g. x[index[ilo]], x[index[ilo+1]]... x[index[ihi]] is ordered
 *
 * Input:
 *   x    : vector of doubles to be ordered from position ilo to ihi
 *   index: vector of integers indexing the values of x
 *   ilo  : first element of x to order
 *   ihi  : last element of x to order
 *   incr : for incr==1 x is returned in increasing order;
 *          incr==-1 in decreasing order
 * Output:
 *   index: rearranged so that x[index[lo]], x[index[lo+1]]...x[index[ihi]]
 *          is ordered
 */
void dindexsort(double *x,
                int *index,
                int ilo,
                int ihi,
                int incr) {
    int pivot;              /* pivot value for partitioning array      */
    int uhi, ulo;           /* indices at ends of unpartitioned region */
    int tempEntry;          /* temporary entry used for swapping       */
    bool sortup, sortlo;    /* indicate if sub-vectors are sorted so
                             * no further subdivision is needed        */

    //assert(x != NULL);
    //assert(index != NULL);

    if (ilo >= ihi) {
        return;
    }

    sortup = sortlo = true;

    /* Select a pivot value */
    pivot = (ilo + ihi) / 2;

    /* Initialize ends of unpartitioned region */
    ulo = ilo;
    uhi = ihi;

    /* While the unpartitioned region is not empty, try to reduce its size */
    while (ulo < uhi) {
        if ((x[index[uhi]] * incr) > (x[index[pivot]] * incr)) {
            /* Check if upper subvector is ordered */
            if ((uhi < ihi) &&
                ((x[index[uhi]] * incr) > (x[index[uhi + 1]] * incr))) {
                sortup = false;
            }

            /* Reduce the size of the unpartitioned region */
            uhi--;
            if ((uhi == pivot) && (ulo < pivot)) {
                tempEntry = index[pivot];
                index[pivot] = index[pivot - 1];
                index[pivot - 1] = tempEntry;
                pivot--;
            }
        } else {
            /* Swap entries at indices ulo and uhi */
            tempEntry = index[ulo];
            index[ulo] = index[uhi];
            index[uhi] = tempEntry;

            if (pivot == ulo) {
                pivot = uhi;
            }

            /* Check if lower subvector is ordered */
            if ((ulo > ilo) &&
                ((x[index[ulo]] * incr) < (x[index[ulo - 1]] * incr))) {
                sortlo = false;
            }

            /* Reduce the size of the unpartitioned region */
            ulo++;
            if ((ulo == pivot) && (uhi > (pivot + 1))) {
                tempEntry = index[pivot];
                index[pivot] = index[pivot + 1];
                index[pivot + 1] = tempEntry;
                pivot++;
            }
        }
    }

    /*
     * Entries from ilo to pivot-1 are < or > pivot and
     * from pivot+1 to ihi are > or < pivot.
     * The two regions can be sorted recursively.
     */
    if ((sortlo == false) && (ilo < (pivot - 1))) {
        dindexsort(x, index, ilo, pivot - 1, incr);
    }
    if ((sortup == false) && (ihi > (pivot + 1))) {
        dindexsort(x, index, pivot + 1, ihi, incr);
    }
}


double runif(void) {
    double x;

    if (cstat_set == 0) {
        setall(is1, is2);
        cstat_set = 1;
    }

    /* assign to double x for conversion */
    x = genunf(0.0, 1.0);
    return (x);
}


/* Returns cdf of normal N(m,s^2) at x */
double pnormC(double y, double m, double s) {
    double cdf, p, mean, sd, bound, x, z;
    /* primitive type conversion */
    x = y;
    mean = m;
    sd = s;
    z = (x - mean) / sd;

    if (z < -20.0) {
        p = 2.753624e-89;
        //p = 2.86e-7F;
    } else if (z > 20.0) {
        p = 1 - 2.753624e-89;
        //p = 0.9999997F;
    } else {
        double q;
        int status;
        int which = 1;
        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);
    }

    cdf = p; /* another primitive type conversion */
    return cdf;
}


/*
 * Density of univariate Normal(m,s^2) evaluated at y.
 * log==1 returns in log-scale.
 */
double dnormC(double y,
              double m,
              double s,
              int logscale) {
    //assert((logscale == 0) || (logscale == 1));

    if (logscale == 1) {
        return (-log(SQ_M_PI_2) - log(s) - 0.5 * (y - m) * (y - m) / (s * s));
    } else {
        return (exp(-0.5 * (y - m) * (y - m) / (s * s)) / (SQ_M_PI_2 * s));
    }
}


/* Joint density of y[0]...y[n-1] under Normal(m,s^2) */
double dnormC_jvec(const double *y,
                   int n,
                   double m,
                   double s,
                   int logscale) {
    int i;
    double ans = 0.0;

    //assert(y != NULL);
    //assert((logscale == 0) || (logscale == 1));

    for (i = 0; i < n; i++) {
        ans += dnormC(y[i], m, s, 1);
    }

    return (logscale == 1) ? ans : exp(ans);
}


/* Returns inv cdf of normal N(m,s^2) at p */
double qnormC(double cdf, double m, double s) {
    double y;

    if ((cdf < 0.0) | (cdf > 1.0)) {
//        errorC("qnormC", "tried inverse cdf with p<0 or p>1", 1);
        /*NOTREACHED*/
    }

    /* par check */
    if (cdf <= 2.753624e-89) {
        y = -20.0 * s + m;
    } else if (cdf >= 0.99999999999999989) {
        y = 8.209536 * s + m;
    } else {
        /* primitive type conversion */
        double p = cdf;
        double q = 1.0 - p;
        double mean = m;
        double sd = s;
        double bound;
        double x;
        int which = 2;
        int status;

        cdfnor(&which, &p, &q, &x, &mean, &sd, &status, &bound);

        y = x; /* another primitive type conversion */
    }

    return y;
}


/*
 * Beta-binomial(alpha,beta) prior probability for a model including
 * k out of p variables.
 */
double bbPrior(int k,
               int p,
               double alpha,
               double beta,
               int logscale) {
    double ans;

    //assert((logscale == 0) || (logscale == 1));

    ans = lnbeta(alpha + k, beta + p - k) - lnbeta(alpha, beta);
    return (logscale == 1) ? ans : exp(ans);
}


/*
 * Complexity prior probability for a model including
 * k out of p variables.
 */
//double complexPrior(int k,
//                    int p,
//                    double priorc,
//                    int logscale) {
//    double priornorm, ans;
//
//    //assert((logscale == 0) || (logscale == 1));
//    priornorm =
//            log(1.0 - 1.0 / pow((double) p, priorc * ((double) p + 1.0))) - log(1.0 - 1.0 / pow((double) p, priorc));
//    ans = lnbeta(1.0 + (double) k, 1.0 + (double) (p - k)) - (priorc * (double) k) * log((double) p) - priornorm;
//
//    return (logscale == 1) ? ans : exp(ans);
//}
//
//
///* Draw from univariate Normal(mu,s^2) */
double rnormC(double mu,
              double s) {
    static bool iset = false;
    static double gset;
    double normdev;

    /* Is a deviate available from a previous invocation? */
    if (iset == false) {
        double fac;
        double rsq;
        double v1;
        double v2;

        do {
            /*
             * Pick two uniform numbers in the square extending
             * from -1 to +1 in each direction
             */
            v1 = 2.0 * runif() - 1.0;
            v2 = 2.0 * runif() - 1.0;
            /* See if they are in the unit circle */
            rsq = (v1 * v1) + (v2 * v2);
        } while (rsq >= 1.0 || rsq == 0.0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        /* Make Box-Muller transformation to get two normal deviates */
        gset = v1 * fac;        /* Save this one for next invocation */
        iset = true;
        normdev = v2 * fac;
    } else {
        iset = false;
        normdev = gset;
    }
    return normdev * s + mu;
}

//Obtain n draws from normal given multiple truncation points, i.e. N(th; m,s) prod_i I(ltrunc[i] < th < rtrunc[i])
//Important note: intervals are assumed to be disjoint and ordered.
// Input
// - n: length of output y
// - ltrunc, rtrunc: vectors [0..ntrunc-1] with start / end of truncation intervals.
// - ntrunc: number of truncatio intervals, i.e. length of vectors ltrunc & rtrunc
// - m: mean of underlying normal
// - s: sd
// Output
// - y: random draws
// - pdfy: joint log-density of y[0..n-1]
void
rnorm_truncMult(double *y, double *pdfy, int *n, double *ltrunc, double *rtrunc, int ntrunc, double *m, double *s) {
    int i, j;
    double u, **p, *cump;
    //Find quantiles for truncation points
    p = dmatrix(0, ntrunc - 1, 0, 1);
    cump = dvector(0, ntrunc);
    cump[0] = 0;
    for (i = 0; i < ntrunc; i++) {
        p[i][0] = pnormC(ltrunc[i], *m, *s);
        p[i][1] = pnormC(rtrunc[i], *m, *s);
        cump[i + 1] = max_xy(cump[i] + 1.0e-30, cump[i] + p[i][1] - p[i][0]);
    }
    //Generate random draws
    (*pdfy) = 0;
    for (i = 0; i < *n; i++) {
        u = runif() * cump[ntrunc];
        j = 0;
        while ((u > cump[j + 1]) && (j < ntrunc - 1)) j++;
        y[i] = qnormC(p[j][0] + u - cump[j], *m, *s);
        (*pdfy) += dnormC(y[i], *m, *s, 1) - log(cump[ntrunc]);
    }
    free_dmatrix(p, 0, ntrunc - 1, 0, 1);
    free_dvector(cump, 0, ntrunc);
}


//Perform 1 Gibbs iteration to sample from z ~ N(alpha, 1) with D[,j] * z[j] restricted outside of interval (lower[j] - D[,-j] %*% z[-j] , upper[j] - D[,-j] %*% z[-j])
//Input:
// - alpha: mean vector [1..p]
// - D: Cholesky decomposition of covariance matrix
// - p: dimensionality of Normal (number of variables)
//Input-Output
// - z: on input contains current draw, on output contains updated draw
// - Dj: on input current value of D %*% z, on output updated according to draw in z
void rtmvnormOutside_Gibbs(double *z, double *Dj, double *alpha, double **D, int p, double *lower, double *upper) {
    int j, k, nrestrict, one = 1;
    double *l, *u, oned = 1;
    l = dvector(1, p);
    u = dvector(1, p);
    for (j = 1; j <= p; j++) {
        //Find truncation points for variable j
        for (k = 1; k <= p; k++) Dj[k] = Dj[k] - D[k][j] * z[j];
        //for (k=1; k<=p; k++) Dj[k]= Dj[k] - D[k][j]*ans[i-1 +(j-1)*n];
        k = 1;
        nrestrict = 0;
        while (k <= p) {
            if (D[k][j] > 0) {
                nrestrict++;
                l[nrestrict] = (lower[k] - Dj[k]) / D[k][j];
                u[nrestrict] = (upper[k] - Dj[k]) / D[k][j];
            } else if (D[k][j] < 0) {
                nrestrict++;
                u[nrestrict] = (lower[k] - Dj[k]) / D[k][j];
                l[nrestrict] = (upper[k] - Dj[k]) / D[k][j];
            }
            k++;
        }
        //Merge excluded regions & define inclusion regions
        if (nrestrict > 0) {
            int *o, jj;
            double *lmod, *umod, lprop;
            o = ivector(1, nrestrict);
            lmod = dvector(1, nrestrict + 1);
            umod = dvector(1, nrestrict + 1);
            for (jj = 1; jj <= nrestrict; jj++) o[jj] = jj;
            dindexsort(l, o, 1, nrestrict, 1); //sort indexes according to values in l
            jj = k = 2;
            lmod[1] = l[o[1]];
            umod[2] = u[o[1]];
            while (jj <= nrestrict) {
                if (u[o[jj]] > umod[k]) {
                    if (l[o[jj]] <= umod[k]) {
                        umod[k] = u[o[jj]];
                    } else {
                        lmod[k] = l[o[jj]];
                        k++;
                        umod[k] = u[o[jj]];
                    }
                }
                jj++;
            }
            umod[1] = -1.0e20;
            lmod[k] = 1.0e20;
            //Draw random variate
            rnorm_truncMult(z + j, &lprop, &one, umod + 1, lmod + 1, k, alpha + j,
                            &oned); //note: umod is lower bound for inclusion intervals (lmod is upper bound)
            for (k = 1; k <= p; k++) Dj[k] += z[j] * D[k][j];
            //rnorm_truncMult(ans+i+(j-1)*n, &lprop, &one, umod+1, lmod+1, k, alpha+j, &oned); //note: umod is lower bound for inclusion intervals (lmod is upper bound)
            //for (k=1; k<=p; k++) Dj[k] += ans[i+(j-1)*n] * D[k][j];
            free_ivector(o, 1, nrestrict);
            free_dvector(lmod, 1, nrestrict + 1);
            free_dvector(umod, 1, nrestrict + 1);
        } else {
            z[j] = rnormC(alpha[j], 1);
        }
    }
    free_dvector(l, 1, p);
    free_dvector(u, 1, p);
}


double dmvtC(const double *y,
             int n,
             const double *mu,
             double **cholsinv,
             double det,
             int nu,
             int logscale) {
    double res = 0.0;
    double t1;
    double t2;
    double normk;

    //assert(y != NULL);
    //assert(mu != NULL);
    //assert(cholsinv != NULL);
    //assert((logscale == 0) || (logscale == 1));

    /* Find (y-mu)' * cholsinv' * cholsinv * (y-mu) */
    {
        int i;
        double *z;
        double *z2;

        z = dvector(1, n);
        z2 = dvector(1, n);
        for (i = 1; i <= n; i++) {
            z[i] = y[i] - mu[i];
        }
        Ax(cholsinv, z, z2, 1, n, 1, n);
        for (i = 1; i <= n; i++) {
            res += z2[i] * z2[i];
        }
        free_dvector(z, 1, n);
        free_dvector(z2, 1, n);
    }

    t2 = 0.5 * nu;
    t1 = t2 + 0.5 * (n + 0.0);
    normk = gamln(&t1) -
            gamln(&t2) -
            0.5 * (n + 0.0) * (log(nu + 0.0) + log(M_PI)) + 0.5 * log(det);

    if (logscale == 1)
        return (normk - t1 * log(1 + res / (nu + 0.0)));
    else
        return (exp(normk) * pow(1.0 + res / (nu + 0.0), -t1));
}

double dmvtC_eigen(const Eigen::VectorXd &y,
                   const Eigen::VectorXd &mu,
                   const Eigen::MatrixXd &cholsinv,
                   double det,
                   int nu,
                   bool logscale) {
    int n = y.size();
    double res;

    // Calculate (y-mu)' * cholsinv' * cholsinv * (y-mu)
    Eigen::VectorXd z = y - mu;
    Eigen::VectorXd z2 = cholsinv * z;
    res = z2.squaredNorm();

    double t2 = 0.5 * nu;
    double t1 = t2 + 0.5 * n;

    // Using our custom gamln function
    double normk = gamln(&t1) -
                   gamln(&t2) -
                   0.5 * n * (std::log(nu) + std::log(M_PI)) +
                   0.5 * std::log(det);

    if (logscale) {
        return normk - t1 * std::log(1 + res / nu);
    } else {
        return std::exp(normk) * std::pow(1.0 + res / nu, -t1);
    }
}

// * Draw from multivar T with n dimensions and nu degrees of freedom
// * Result is stored in y[1..n].
// *     mu is the location parameter
// *     chols is the Cholesky decomposition of the covariance matrix.
// * That is, the covariance is s*nu/(nu-2) and s=chols*chols'
// * and nu are the degrees of freedom
// * Note: both y and mu should have length n, and s should be an n*n matrix.
// * The routine doesn't check it.
// *
// * Example:
// *   choldc(s,n,chols,posdef); //compute cholesky decomposition
// *   rmvtC(y,n,mu,chols,nu); //generate random variate
// */
void rmvtC(double *y,
           int n,
           const double *mu,
           double **chols,
           int nu) {
    int i;
    double x;
    double *z;

    //assert(y != NULL);
    //assert(mu != NULL);
    //assert(chols != NULL);

    /* Draw from chi-square with nu degrees of freedom */
    x = sqrt(nu / gengam(0.5, nu / 2.0));

    /* Multiple n indep normal draws by the common chi-square */
    z = dvector(1, n);
    for (i = 1; i <= n; i++) {
        z[i] = x * rnormC(0, 1);
    }
    /* Compute mu + chols*z */
    Ax_plus_y(chols, z, mu, y, 1, n);

    free_dvector(z, 1, n);
}

void rmvtC_eigen(Eigen::VectorXd &y, const Eigen::VectorXd &mu, const Eigen::MatrixXd &chols, int nu) {
    int n = mu.size();
    double x = std::sqrt(nu / gengam(0.5, nu / 2.0));
    Eigen::VectorXd z = Eigen::VectorXd::NullaryExpr(n, [&]() { return x * rnormC(0, 1); });
    y = mu + chols * z;
}

//Univariate iMOM prior
double dimom(double y, double m, double tau, double phi, int logscale) {
    double y2, ans;
    y2 = (y - m) * (y - m);
    ans = .5 * (log(tau) + log(phi)) - .5 * LOG_M_PI - log(y2) - tau * phi / y2;
    if (logscale == 0) ans = exp(ans);
    return (ans);
}


void
rnlpPost_lm(double *ans, int niter, int burnin, int thinning, double *y, int *x, int n, int p, int r, double tau,
            double a_phi, double b_phi, int prior, double *thini, double phiini) {
    bool posdef;
    int i, j, k, isave, nsave;
    double *m, *mortho, *alpha, **S, **Sinv, **cholSinv, **inv_cholSinv, **K, **D, tauinv =
            1.0 / tau, *Xty, *thcur, phicur, phinew, sqrtphi, th2sum, th2invsum, apost, bpost, *linpred, ssr;
    //Pre-compute stuff
    nsave = (int) floor((niter - burnin + .0) / (thinning + .0));
    m = dvector(1, p);
    mortho = dvector(1, p);
    alpha = dvector(1, p);
    thcur = dvector(1, p);
    linpred = dvector(0, n - 1);
    S = dmatrix(1, p, 1, p);
    Sinv = dmatrix(1, p, 1, p);
    cholSinv = dmatrix(1, p, 1, p);
    inv_cholSinv = dmatrix(1, p, 1, p);
    K = dmatrix(1, p, 1, p);
    D = dmatrix(1, p, 1, p);

    AvectBvec_int(x, n, p, x, n, p, S); //S= t(x) %*% x + 1/tau
    for (i = 1; i <= p; i++) {
        S[i][i] += tauinv;
    }

//    printMatrix("S", S, p, p);

    inv_posdef(S, p, Sinv, &posdef);
//    printMatrix("Sinv", Sinv, p, p);

    choldc(Sinv, p, cholSinv, &posdef);
//    printMatrix("cholSinv", cholSinv, p, p);

    choldc_inv(Sinv, p, inv_cholSinv, &posdef); //inverse of chol(Sinv)

    Xty = dvector(1, p);
    Atvecx_int(x, y, Xty + 1, 0, p - 1, 0, n - 1); //m= solve(S) %*% t(x) %*% y
    Ax(Sinv, Xty, m, 1, p, 1, p);
    Ax(inv_cholSinv, m, mortho, 1, p, 1, p);
    free_dvector(Xty, 1, p);

    if (prior == 0) apost = .5 * (a_phi + n + 3 * p);
    else if (prior == 1)
        apost = .5 * (a_phi + n - p);
    else apost = .5 * (a_phi + n + p);
    //Initialize
    th2sum = 0;


    phicur = phiini;
    sqrtphi = sqrt(phiini);
//    phicur = sqrtphi = 1.0;
    for (j = 1; j <= p; j++) {
        thcur[j] = thini[j - 1];
//        std::cout << thcur[j] << " ";
//        thcur[j] = m[j];
        th2sum += thcur[j] * thcur[j];
    }

//    std::cout << std::endl;


//    double a = thcur[1];

    //Ax(cholSinv, thcur, Dthcur, 1, p, 1, p);
    isave = 0;
    // TODO loop wahrscheinlich raus
    for (i = 1; i <= niter; i++) {
        //for (j=1; j<=p; j++) Dthcur[j] = Dthcur[j] / sqrtphi;
        Avecx_int(x, thcur + 1, linpred, 0, n - 1, 0, p - 1);
        ssr = 0;
        for (j = 0; j < n; j++) {
            ssr += pow(y[j] - linpred[j], 2.0);
        }
//        if (prior == 0) {
//            bpost = .5 * (b_phi + th2sum / tau + ssr);
//            phicur = 1.0 / rgammaC(apost, bpost);
//            sqrtphi = sqrt(phicur);
//        } else {
//            if (prior == 1) bpost = .5 * (b_phi + ssr); else bpost = .5 * (b_phi + th2sum / tau + ssr);
//            phinew = 1.0 / rgammaC(apost, bpost);
//            th2invsum = 0;
//            for (j = 1; j <= p; j++) th2invsum += 1 / (thcur[j] * thcur[j]);
//            if (runif() < exp((phicur - phinew) * tau * th2invsum)) {
//                phicur = phinew;
//                sqrtphi = sqrt(phicur);
//            }
//        }
        for (j = 1; j <= p; j++) {
            alpha[j] = mortho[j] / sqrtphi;
            //Dthcur[j]= Dthcur[j] * sqrtphi;
            for (k = 1; k <= j; k++) {
                D[j][k] = cholSinv[j][k] * sqrtphi;
                K[j][k] = inv_cholSinv[j][k] / sqrtphi;
            }
        }

        // Print the matrices
//        printMatrix("D matrix", D, p, p);
//        printMatrix("K matrix", K, p, p);
//

        rnlp_Gibbs(thcur, p, alpha, D, K, &tau, &phicur, r, prior);

//        for (j = 1; j <= p; j++) {
//            std::cout << thcur[j] << " ";
//        }
//
//        std::cout << std::endl;


        if (i > burnin && ((i - burnin) % thinning) == 0) {
            for (j = 1; j <= p; j++) {
                ans[j - 1] = thcur[j];
                if (isnan(thcur[j])) {
                    std::cout << "";
                }
            }
//            ans[isave + p * nsave] = phicur;
//            isave++; // Muss erst checken ob das rein oder raus gehört.
        }
    }
    free_dvector(m, 1, p);
    free_dvector(mortho, 1, p);
    free_dvector(alpha, 1, p);
    free_dvector(thcur, 1, p);
    free_dvector(linpred, 0, n - 1);
    free_dmatrix(S, 1, p, 1, p);
    free_dmatrix(Sinv, 1, p, 1, p);
    free_dmatrix(cholSinv, 1, p, 1, p);
    free_dmatrix(inv_cholSinv, 1, p, 1, p);
    free_dmatrix(K, 1, p, 1, p);
    free_dmatrix(D, 1, p, 1, p);
}

//Single Gibbs update (th,l) ~ N(th;m,S) * prod g(th[i]) >l[i]
//Input
// - p: dimensionality of th (number of variables)
// - m: m[1..p] is mean of Normal kernel
// - cholS: cholS[1..p][1..p] is Cholesky decomp of covariance
// - K: inverse of cholS
// - tau: value of tau (prior dispersion)
// - phi: value of phi (residual variance)
// - r: power parameter is 2*r
// - prior: prior==0 for MOM, prior==1 for iMOM, prior==2 for eMOM
//Input-Output
// - th: at input th[1..p] is current value of th; at output the updated value
void rnlp_Gibbs(double *th, int p, double *m, double **cholS, double **K, double *tau, double *phi, int r, int prior) {
    int i;
    double *lower, *upper, *l, *z, upperb;
    lower = dvector(1, p);
    upper = dvector(1, p);
    l = dvector(1, p);
    z = dvector(1, p);
    //Update l
    if (prior == 0) {
        for (i = 1; i <= p; i++) {
            upperb = pen_mom(th + i, phi, tau, r);
            l[i] = runif() * upperb;  //l[i]= runif() * pow(th[i]*th[i] / (phi*tau), r + .0);
            if (r == 1) { upper[i] = sqrt(l[i] * (*tau) * (*phi)); }
            else {
                upper[i] = pow(l[i] * (*tau) * (*phi), 1.0 / (2.0 * r));
            }
            lower[i] = -upper[i];
        }
    } else if (prior == 1) {
        for (i = 1; i <= p; i++) {
            upperb = pen_imom(th + i, phi, tau, 1);
            l[i] = log(runif()) + upperb;
            upper[i] = invpen_imom_sandwich(l + i, phi, tau);
            lower[i] = -upper[i];
        }
    } else if (prior == 2) {
        for (i = 1; i <= p; i++) {
            upperb = pen_emom(th + i, phi, tau, 1);
            l[i] = runif() * exp(upperb);  //l[i]= runif() * exp(sqrt(2) + tau*phi/th[i]^2)
            upper[i] = sqrt(fabs((*tau) * (*phi) / (log(l[i]) - sqrt(2.0))));
            lower[i] = -upper[i];
        }
    }
    //Update th, cholSth
    Ax(K, th, z, 1, p, 1, p); //z= K th;
    rtmvnormOutside_Gibbs(z, th, m, cholS, p, lower, upper);
    Ax(cholS, z, th, 1, p, 1, p); //th= D z
    free_dvector(lower, 1, p);
    free_dvector(upper, 1, p);
    free_dvector(l, 1, p);
    free_dvector(z, 1, p);
}


////Evaluates MOM prior penalty, i.e. (th^2 / (phi*tau))^r
double pen_mom(double *th, double *phi, double *tau, int r) {
    double ans;
    ans = pow(th[0] * th[0] / ((*phi) * (*tau)), r + .0);
    return ans;
}

//
////Evaluates eMOM prior penalty, i.e. exp(-sqrt(2)*tau*phi/th^2)
double pen_emom(double *th, double *phi, double *tau, int logscale) {
    double ans;
    ans = sqrt(2.0) - (*tau) * (*phi) / (th[0] * th[0]);
    if (logscale == 0) ans = exp(ans);
    return ans;
}

//Evaluates iMOM prior penalty, i.e. dimom(th) / dnorm(th)
double pen_imom(double *th, double *phi, double *tau, int logscale) {
    double ans;
    ans = dimom(*th, 0, *tau, *phi, 1) - dnormC(*th, 0, sqrt((*tau) * (*phi)), 1);
    if (logscale == 0) ans = exp(ans);
    return (ans);
}


// Uses an initial guess to bound the solution. Then uses a sandwhich approach based on recursive linear interpolation
double invpen_imom_sandwich(double *loglambda, double *phi, double *tau) {
    int i, maxiter = 50;
    double b, d, zcur, thcur, fcur, zlow, thlow, flow, zup, thup, fup, err, ftol = 1.0e-5;
    //Initial guess
    b = .5 * (log((*tau) * (*tau)) + 2.0 * log(*phi) + log(2.0)) - (*loglambda);
    d = sqrt(b * b + 2.0);
    zcur = (*tau) * (*phi) * (-b + d);
    thcur = sqrt(zcur);
    fcur = pen_imom(&thcur, phi, tau, 1);  //log-penalty at zcur
    //Lower & upper bound
    if (fcur >= (*loglambda)) {
        zlow = .8 * .8 * zcur;
        thlow = sqrt(zlow);
        flow = pen_imom(&thlow, phi, tau, 1);
        while (flow >= (*loglambda)) {
            zcur = zlow;
            thcur = thlow;
            fcur = flow;
            zlow = .8 * .8 * zcur;
            thlow = sqrt(zlow);
            flow = pen_imom(&thlow, phi, tau, 1);
        }
        zup = zcur;
        thup = thcur;
        fup = fcur;
    } else {
        zup = 1.2 * 1.2 * zcur;
        thup = sqrt(zup);
        fup = pen_imom(&thup, phi, tau, 1);
        while (fup <= (*loglambda)) {
            zcur = zup;
            thcur = thup;
            fcur = fup;
            zup = 1.2 * 1.2 * zcur;
            thup = sqrt(zup);
            fup = pen_imom(&thup, phi, tau, 1);
        }
        zlow = zcur;
        thlow = thcur;
        flow = fcur;
    }
    //Search by sandwich linear interpolation
    err = fcur - *loglambda;
    i = 1;
    while ((i < maxiter) && (fabs(err) > ftol)) {
        b = (fup - flow) / (zup - zlow);
        zcur = zlow + ((*loglambda) - flow) / b; //approx is flow + b*(z-zlow)= loglambda
        thcur = sqrt(zcur);
        fcur = pen_imom(&thcur, phi, tau, 1);
        err = fcur - *loglambda;
        if (err > 0) {
            zup = zcur;
            fup = fcur;
        } else {
            zlow = zcur;
            flow = fcur;
        }
        i++;
    }
    return thcur;
}


/************************************************************************
                       MORE RANDOM VARIATE STUFF
************************************************************************/

/*
 * Truncates a double precision number to an integer.
 *     a - number to be truncated
 */
double fifdint(double a) {
    long temp = (long) (a);
    return (double) (temp);
}


/**********************************************************************

      void cdfnor(int *which,double *p,double *q,double *x,double *mean,
            double *sd,int *status,double *bound)

               Cumulative Distribution Function
               NORmal distribution


                              Function


     Calculates any one parameter of the normal
     distribution given values for the others.


                              Arguments


     WHICH  --> Integer indicating  which of the  next  parameter
     values is to be calculated using values  of the others.
     Legal range: 1..4
               iwhich = 1 : Calculate P and Q from X,MEAN and SD
               iwhich = 2 : Calculate X from P,Q,MEAN and SD
               iwhich = 3 : Calculate MEAN from P,Q,X and SD
               iwhich = 4 : Calculate SD from P,Q,X and MEAN

     P <--> The integral from -infinity to X of the normal density.
            Input range: (0,1].

     Q <--> 1-P.
            Input range: (0, 1].
            P + Q = 1.0.

     X < --> Upper limit of integration of the normal-density.
             Input range: ( -infinity, +infinity)

     MEAN <--> The mean of the normal density.
               Input range: (-infinity, +infinity)

     SD <--> Standard Deviation of the normal density.
             Input range: (0, +infinity).

     STATUS <-- 0 if calculation completed correctly
               -I if input parameter number I is out of range
                1 if answer appears to be lower than lowest
                  search bound
                2 if answer appears to be higher than greatest
                  search bound
                3 if P + Q .ne. 1

     BOUND <-- Undefined if STATUS is 0

               Bound exceeded by parameter number I if STATUS
               is negative.

               Lower search bound if STATUS is 1.

               Upper search bound if STATUS is 2.


                              Method




     A slightly modified version of ANORM from

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     is used to calulate the  cumulative standard normal distribution.

     The rational functions from pages  90-95  of Kennedy and Gentle,
     Statistical  Computing,  Marcel  Dekker, NY,  1980 are  used  as
     starting values to Newton's Iterations which compute the inverse
     standard normal.  Therefore no  searches  are necessary for  any
     parameter.

     For X < -15, the asymptotic expansion for the normal is used  as
     the starting value in finding the inverse standard normal.
     This is formula 26.2.12 of Abramowitz and Stegun.


                              Note


      The normal density is proportional to
      exp( - 0.5 * (( X - MEAN)/SD)**2)

**********************************************************************/
void cdfnor(int *which,
            double *p,
            double *q,
            double *x,
            double *mean,
            double *sd,
            int *status,
            double *bound) {
    static int K1 = 1;
    static double z, pq;
/*
     ..
     .. Executable Statements ..
*/
/*
     Check arguments
*/
    *status = 0;
    if (!(*which < 1 || *which > 4)) goto S30;
    if (!(*which < 1)) goto S10;
    *bound = 1.0e0;
    goto S20;
    S10:
    *bound = 4.0e0;
    S20:
    *status = -1;
    return;
    S30:
    if (*which == 1) goto S70;
/*
     P
*/
    if (!(*p <= 0.0e0 || *p > 1.0e0)) goto S60;
    if (!(*p <= 0.0e0)) goto S40;
    *bound = 0.0e0;
    goto S50;
    S40:
    *bound = 1.0e0;
    S50:
    *status = -2;
    return;
    S70:
    S60:
    if (*which == 1) goto S110;
/*
     Q
*/
    if (!(*q <= 0.0e0 || *q > 1.0e0)) goto S100;
    if (!(*q <= 0.0e0)) goto S80;
    *bound = 0.0e0;
    goto S90;
    S80:
    *bound = 1.0e0;
    S90:
    *status = -3;
    return;
    S110:
    S100:
    if (*which == 1) goto S150;
/*
     P + Q
*/
    pq = *p + *q;
    if (!(fabs(pq - 0.5e0 - 0.5e0) > 3.0e0 * spmpar(&K1))) goto S140;
    if (!(pq < 0.0e0)) goto S120;
    *bound = 0.0e0;
    goto S130;
    S120:
    *bound = 1.0e0;
    S130:
    *status = 3;
    return;
    S150:
    S140:
    if (*which == 4) goto S170;
/*
     SD
*/
    if (!(*sd <= 0.0e0)) goto S160;
    *bound = 0.0e0;
    *status = -6;
    return;
    S170:
    S160:
/*
     Calculate ANSWERS
*/
    if (1 == *which) {
/*
     Computing P
*/
        z = (*x - *mean) / *sd;
        cumnor(&z, p, q);
    } else if (2 == *which) {
/*
     Computing X
*/
        z = dinvnr(p, q);
        *x = *sd * z + *mean;
    } else if (3 == *which) {
/*
     Computing the MEAN
*/
        z = dinvnr(p, q);
        *mean = *x - *sd * z;
    } else if (4 == *which) {
/*
     Computing SD
*/
        z = dinvnr(p, q);
        *sd = (*x - *mean) / z;
    }
    return;
}


/*
-----------------------------------------------------------------------

     SPMPAR PROVIDES THE SINGLE PRECISION MACHINE CONSTANTS FOR
     THE COMPUTER BEING USED. IT IS ASSUMED THAT THE ARGUMENT
     I IS AN INTEGER HAVING ONE OF THE VALUES 1, 2, OR 3. IF THE
     SINGLE PRECISION ARITHMETIC BEING USED HAS M BASE B DIGITS AND
     ITS SMALLEST AND LARGEST EXPONENTS ARE EMIN AND EMAX, THEN

        SPMPAR(1) = B**(1 - M), THE MACHINE PRECISION,

        SPMPAR(2) = B**(EMIN - 1), THE SMALLEST MAGNITUDE,

        SPMPAR(3) = B**EMAX*(1 - B**(-M)), THE LARGEST MAGNITUDE.

-----------------------------------------------------------------------
     WRITTEN BY
        ALFRED H. MORRIS, JR.
        NAVAL SURFACE WARFARE CENTER
        DAHLGREN VIRGINIA
-----------------------------------------------------------------------
-----------------------------------------------------------------------
     MODIFIED BY BARRY W. BROWN TO RETURN DOUBLE PRECISION MACHINE
     CONSTANTS FOR THE COMPUTER BEING USED.  THIS MODIFICATION WAS
     MADE AS PART OF CONVERTING BRATIO TO DOUBLE PRECISION
-----------------------------------------------------------------------
*/
double spmpar(int *i) {
    static int K1 = 4;
    static int K2 = 8;
    static int K3 = 9;
    static int K4 = 10;
    static double spmpar, b, binv, bm1, one, w, z;
    static int emax, emin, ibeta, m;
/*
     ..
     .. Executable Statements ..
*/
    if (*i > 1) goto S10;
    b = ipmpar(&K1);
    m = ipmpar(&K2);
    spmpar = pow(b, (double) (1 - m));
    return spmpar;
    S10:
    if (*i > 2) goto S20;
    b = ipmpar(&K1);
    emin = ipmpar(&K3);
    one = 1.0;
    binv = one / b;
    w = pow(b, (double) (emin + 2));
    spmpar = w * binv * binv * binv;
    return spmpar;
    S20:
    ibeta = ipmpar(&K1);
    m = ipmpar(&K2);
    emax = ipmpar(&K4);
    b = ibeta;
    bm1 = ibeta - 1;
    one = 1.0;
    z = pow(b, (double) (m - 1));
    w = ((z - one) * b + bm1) / (b * z);
    z = pow(b, (double) (emax - 2));
    spmpar = w * z * b * b;
    return spmpar;
}


/*
**********************************************************************

     void cumnor(double *arg,double *result,double *ccum)


                              Function


     Computes the cumulative  of    the  normal   distribution,   i.e.,
     the integral from -infinity to x of
          (1/sqrt(2*pi)) exp(-u*u/2) du

     X --> Upper limit of integration.
                                        X is DOUBLE PRECISION

     RESULT <-- Cumulative normal distribution.
                                        RESULT is DOUBLE PRECISION

     CCUM <-- Compliment of Cumulative normal distribution.
                                        CCUM is DOUBLE PRECISION

     Renaming of function ANORM from:

     Cody, W.D. (1993). "ALGORITHM 715: SPECFUN - A Portabel FORTRAN
     Package of Special Function Routines and Test Drivers"
     acm Transactions on Mathematical Software. 19, 22-32.

     with slight modifications to return ccum and to deal with
     machine constants.

**********************************************************************
  Original Comments:
------------------------------------------------------------------

 This function evaluates the normal distribution function:

                              / x
                     1       |       -t*t/2
          P(x) = ----------- |      e       dt
                 sqrt(2 pi)  |
                             /-oo

   The main computation evaluates near-minimax approximations
   derived from those in "Rational Chebyshev approximations for
   the error function" by W. J. Cody, Math. Comp., 1969, 631-637.
   This transportable program uses rational functions that
   theoretically approximate the normal distribution function to
   at least 18 significant decimal digits.  The accuracy achieved
   depends on the arithmetic system, the compiler, the intrinsic
   functions, and proper selection of the machine-dependent
   constants.

*******************************************************************
*******************************************************************

 Explanation of machine-dependent constants.

   MIN   = smallest machine representable number.

   EPS   = argument below which anorm(x) may be represented by
           0.5  and above which  x*x  will not underflow.
           A conservative value is the largest machine number X
           such that   1.0 + X = 1.0   to machine precision.
*******************************************************************
*******************************************************************

 Error returns

  The program returns  ANORM = 0     for  ARG .LE. XLOW.


 Intrinsic functions required are:

     ABS, AINT, EXP


  Author: W. J. Cody
          Mathematics and Computer Science Division
          Argonne National Laboratory
          Argonne, IL 60439

  Latest modification: March 15, 1992

------------------------------------------------------------------
*/
void cumnor(double *arg,
            double *result,
            double *ccum) {
    static double a[5] = {
            2.2352520354606839287e00, 1.6102823106855587881e02, 1.0676894854603709582e03,
            1.8154981253343561249e04, 6.5682337918207449113e-2
    };
    static double b[4] = {
            4.7202581904688241870e01, 9.7609855173777669322e02, 1.0260932208618978205e04,
            4.5507789335026729956e04
    };
    static double c[9] = {
            3.9894151208813466764e-1, 8.8831497943883759412e00, 9.3506656132177855979e01,
            5.9727027639480026226e02, 2.4945375852903726711e03, 6.8481904505362823326e03,
            1.1602651437647350124e04, 9.8427148383839780218e03, 1.0765576773720192317e-8
    };
    static double d[8] = {
            2.2266688044328115691e01, 2.3538790178262499861e02, 1.5193775994075548050e03,
            6.4855582982667607550e03, 1.8615571640885098091e04, 3.4900952721145977266e04,
            3.8912003286093271411e04, 1.9685429676859990727e04
    };
    static double half = 0.5e0;
    static double p[6] = {
            2.1589853405795699e-1, 1.274011611602473639e-1, 2.2235277870649807e-2,
            1.421619193227893466e-3, 2.9112874951168792e-5, 2.307344176494017303e-2
    };
    static double one = 1.0e0;
    static double q[5] = {
            1.28426009614491121e00, 4.68238212480865118e-1, 6.59881378689285515e-2,
            3.78239633202758244e-3, 7.29751555083966205e-5
    };
    static double sixten = 1.60e0;
    static double sqrpi = 3.9894228040143267794e-1;
    static double thrsh = 0.66291e0;
    static double root32 = 5.656854248e0;
    static double zero = 0.0e0;
    static int K1 = 1;
    static int K2 = 2;
    static int i;
    static double del, eps, temp, x, xden, xnum, y, xsq, min;
/*
------------------------------------------------------------------
  Machine dependent constants
------------------------------------------------------------------
*/
    eps = spmpar(&K1) * 0.5e0;
    min = spmpar(&K2);
    x = *arg;
    y = fabs(x);
    if (y <= thrsh) {
/*
------------------------------------------------------------------
  Evaluate  anorm  for  |X| <= 0.66291
------------------------------------------------------------------
*/
        xsq = zero;
        if (y > eps) xsq = x * x;
        xnum = a[4] * xsq;
        xden = xsq;
        for (i = 0; i < 3; i++) {
            xnum = (xnum + a[i]) * xsq;
            xden = (xden + b[i]) * xsq;
        }
        *result = x * (xnum + a[3]) / (xden + b[3]);
        temp = *result;
        *result = half + temp;
        *ccum = half - temp;
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for 0.66291 <= |X| <= sqrt(32)
------------------------------------------------------------------
*/
    else if (y <= root32) {
        xnum = c[8] * y;
        xden = y;
        for (i = 0; i < 7; i++) {
            xnum = (xnum + c[i]) * y;
            xden = (xden + d[i]) * y;
        }
        *result = (xnum + c[7]) / (xden + d[7]);
        xsq = fifdint(y * sixten) / sixten;
        del = (y - xsq) * (y + xsq);
        *result = exp(-(xsq * xsq * half)) * exp(-(del * half)) * *result;
        *ccum = one - *result;
        if (x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
/*
------------------------------------------------------------------
  Evaluate  anorm  for |X| > sqrt(32)
------------------------------------------------------------------
*/
    else {
        *result = zero;
        xsq = one / (x * x);
        xnum = p[5] * xsq;
        xden = xsq;
        for (i = 0; i < 4; i++) {
            xnum = (xnum + p[i]) * xsq;
            xden = (xden + q[i]) * xsq;
        }
        *result = xsq * (xnum + p[4]) / (xden + q[4]);
        *result = (sqrpi - *result) / y;
        xsq = fifdint(x * sixten) / sixten;
        del = (x - xsq) * (x + xsq);
        *result = exp(-(xsq * xsq * half)) * exp(-(del * half)) * *result;
        *ccum = one - *result;
        if (x > zero) {
            temp = *result;
            *result = *ccum;
            *ccum = temp;
        }
    }
    if (*result < min) *result = 0.0e0;
/*
------------------------------------------------------------------
  Fix up for negative argument, erf, etc.
------------------------------------------------------------------
----------Last card of ANORM ----------
*/
    if (*ccum < min) *ccum = 0.0e0;
}


/*
**********************************************************************

     double dinvnr(double *p,double *q)
     Double precision NoRmal distribution INVerse


                              Function


     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


                              Arguments


     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION

     Q --> 1-P
                    P is DOUBLE PRECISION


                              Method


     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980 is used as a start
     value for the Newton method of finding roots.


                              Note


     If P or Q .lt. machine EPS returns +/- DINVNR(EPS)

**********************************************************************
*/
double dinvnr(double *p,
              double *q) {
#define maxit 100
#define eps 1.0e-13
#define r2pi 0.3989422804014326e0
#define nhalf -0.5e0
#define dennor(x) (r2pi*exp(nhalf*(x)*(x)))
    static double dinvnr, strtx, xcur, cum, ccum, pp, dx;
    static int i;
    static unsigned long qporq;
/*
     ..
     .. Executable Statements ..
*/
/*
     FIND MINIMUM OF P AND Q
*/
    qporq = *p <= *q;
    if (!qporq) goto S10;
    pp = *p;
    goto S20;
    S10:
    pp = *q;
    S20:
/*
     INITIALIZATION STEP
*/
    strtx = stvaln(&pp);
    xcur = strtx;
/*
     NEWTON INTERATIONS
*/
    for (i = 1; i <= maxit; i++) {
        cumnor(&xcur, &cum, &ccum);
        dx = (cum - pp) / dennor(xcur);
        xcur -= dx;
        if (fabs(dx / xcur) < eps) goto S40;
    }
    dinvnr = strtx;
/*
     IF WE GET HERE, NEWTON HAS FAILED
*/
    if (!qporq) dinvnr = -dinvnr;
    return dinvnr;
    S40:
/*
     IF WE GET HERE, NEWTON HAS SUCCEDED
*/
    dinvnr = xcur;
    if (!qporq) dinvnr = -dinvnr;
    return dinvnr;
#undef maxit
#undef eps
#undef r2pi
#undef nhalf
#undef dennor
}


/*
**********************************************************************

     double stvaln(double *p)
                    STarting VALue for Neton-Raphon
                calculation of Normal distribution Inverse


                              Function


     Returns X  such that CUMNOR(X)  =   P,  i.e., the  integral from -
     infinity to X of (1/SQRT(2*PI)) EXP(-U*U/2) dU is P


                              Arguments


     P --> The probability whose normal deviate is sought.
                    P is DOUBLE PRECISION


                              Method


     The  rational   function   on  page 95    of Kennedy  and  Gentle,
     Statistical Computing, Marcel Dekker, NY , 1980.

**********************************************************************
*/
double stvaln(double *p) {
    static double xden[5] = {
            0.993484626060e-1, 0.588581570495e0, 0.531103462366e0, 0.103537752850e0,
            0.38560700634e-2
    };
    static double xnum[5] = {
            -0.322232431088e0, -1.000000000000e0, -0.342242088547e0, -0.204231210245e-1,
            -0.453642210148e-4
    };
    static int K1 = 5;
    static double stvaln, sign, y, z;
/*
     ..
     .. Executable Statements ..
*/
    if (!(*p <= 0.5e0)) goto S10;
    sign = -1.0e0;
    z = *p;
    goto S20;
    S10:
    sign = 1.0e0;
    z = 1.0e0 - *p;
    S20:
    y = sqrt(-(2.0e0 * log(z)));
    stvaln = y + devlpl(xnum, &K1, &y) / devlpl(xden, &K1, &y);
    stvaln = sign * stvaln;
    return stvaln;
}


/*
**********************************************************************

     double devlpl(double a[],int *n,double *x)
              Double precision EVALuate a PoLynomial at X


                              Function


     returns
          A(1) + A(2)*X + ... + A(N)*X**(N-1)


                              Arguments


     A --> Array of coefficients of the polynomial.
                                        A is DOUBLE PRECISION(N)

     N --> Length of A, also degree of polynomial - 1.
                                        N is INTEGER

     X --> Point at which the polynomial is to be evaluated.
                                        X is DOUBLE PRECISION

**********************************************************************
*/
double devlpl(double a[],
              int *n,
              double *x) {
    static double devlpl, term;
    static int i;
/*
     ..
     .. Executable Statements ..
*/
    term = a[*n - 1];
    for (i = *n - 1 - 1; i >= 0; i--) term = a[i] + term * *x;
    devlpl = term;
    return devlpl;
}


/*
-----------------------------------------------------------------------

     IPMPAR PROVIDES THE INTEGER MACHINE CONSTANTS FOR THE COMPUTER
     THAT IS USED. IT IS ASSUMED THAT THE ARGUMENT I IS AN INTEGER
     HAVING ONE OF THE VALUES 1-10. IPMPAR(I) HAS THE VALUE ...

  INTEGERS.

     ASSUME INTEGERS ARE REPRESENTED IN THE N-DIGIT, BASE-A FORM

               SIGN ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )

               WHERE 0 .LE. X(I) .LT. A FOR I=0,...,N-1.

     IPMPAR(1) = A, THE BASE.

     IPMPAR(2) = N, THE NUMBER OF BASE-A DIGITS.

     IPMPAR(3) = A**N - 1, THE LARGEST MAGNITUDE.

  FLOATING-POINT NUMBERS.

     IT IS ASSUMED THAT THE SINGLE AND DOUBLE PRECISION FLOATING
     POINT ARITHMETICS HAVE THE SAME BASE, SAY B, AND THAT THE
     NONZERO NUMBERS ARE REPRESENTED IN THE FORM

               SIGN (B**E) * (X(1)/B + ... + X(M)/B**M)

               WHERE X(I) = 0,1,...,B-1 FOR I=1,...,M,
               X(1) .GE. 1, AND EMIN .LE. E .LE. EMAX.

     IPMPAR(4) = B, THE BASE.

  SINGLE-PRECISION

     IPMPAR(5) = M, THE NUMBER OF BASE-B DIGITS.

     IPMPAR(6) = EMIN, THE SMALLEST EXPONENT E.

     IPMPAR(7) = EMAX, THE LARGEST EXPONENT E.

  DOUBLE-PRECISION

     IPMPAR(8) = M, THE NUMBER OF BASE-B DIGITS.

     IPMPAR(9) = EMIN, THE SMALLEST EXPONENT E.

     IPMPAR(10) = EMAX, THE LARGEST EXPONENT E.

-----------------------------------------------------------------------

     TO DEFINE THIS FUNCTION FOR THE COMPUTER BEING USED REMOVE
     THE COMMENT DELIMITORS FROM THE DEFINITIONS DIRECTLY BELOW THE NAME
     OF THE MACHINE

-----------------------------------------------------------------------

     IPMPAR IS AN ADAPTATION OF THE FUNCTION I1MACH, WRITTEN BY
     P.A. FOX, A.D. HALL, AND N.L. SCHRYER (BELL LABORATORIES).
     IPMPAR WAS FORMED BY A.H. MORRIS (NSWC). THE CONSTANTS ARE
     FROM BELL LABORATORIES, NSWC, AND OTHER SOURCES.

-----------------------------------------------------------------------
     .. Scalar Arguments ..
*/
int ipmpar(int *i) {
    static int imach[11];
    static int ipmpar;
/*     MACHINE CONSTANTS FOR AMDAHL MACHINES. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
       PC 7300, AND AT&T 6300. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 33;
   imach[3] = 8589934591;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -256;
   imach[7] = 255;
   imach[8] = 60;
   imach[9] = -256;
   imach[10] = 255;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -50;
   imach[10] = 76;
*/
/*     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS. */
/*
   imach[1] = 2;
   imach[2] = 39;
   imach[3] = 549755813887;
   imach[4] = 8;
   imach[5] = 13;
   imach[6] = -50;
   imach[7] = 76;
   imach[8] = 26;
   imach[9] = -32754;
   imach[10] = 32780;
*/
/*     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 48;
   imach[3] = 281474976710655;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -974;
   imach[7] = 1070;
   imach[8] = 95;
   imach[9] = -926;
   imach[10] = 1070;
*/
/*     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
       ARITHMETIC (NOS/VE OPERATING SYSTEM). */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 48;
   imach[6] = -4096;
   imach[7] = 4095;
   imach[8] = 96;
   imach[9] = -4096;
   imach[10] = 4095;
*/
/*     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3. */
/*
   imach[1] = 2;
   imach[2] = 63;
   imach[3] = 9223372036854775807;
   imach[4] = 2;
   imach[5] = 47;
   imach[6] = -8189;
   imach[7] = 8190;
   imach[8] = 94;
   imach[9] = -8099;
   imach[10] = 8190;
*/
/*     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200. */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE HARRIS 220. */
/*
   imach[1] = 2;
   imach[2] = 23;
   imach[3] = 8388607;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 38;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HONEYWELL 600/6000
       AND DPS 8/70 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 63;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       3 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 39;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 2100
       4 WORD DOUBLE PRECISION OPTION WITH FTN4 */
/*
   imach[1] = 2;
   imach[2] = 15;
   imach[3] = 32767;
   imach[4] = 2;
   imach[5] = 23;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 55;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE HP 9000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -126;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
       5/7/9 AND THE SEL SYSTEMS 85/86. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 16;
   imach[5] = 6;
   imach[6] = -64;
   imach[7] = 63;
   imach[8] = 14;
   imach[9] = -64;
   imach[10] = 63;
*/
/*     MACHINE CONSTANTS FOR THE IBM PC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
       MACFORTRAN II. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 54;
   imach[9] = -101;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 62;
   imach[9] = -128;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
       32-BIT INTEGER ARITHMETIC. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
/*     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
       SERIES (MIPS R3000 PROCESSOR). */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -125;
   imach[7] = 128;
   imach[8] = 53;
   imach[9] = -1021;
   imach[10] = 1024;
*/
/*     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300). */

    imach[1] = 2;
    imach[2] = 31;
    imach[3] = 2147483647;
    imach[4] = 2;
    imach[5] = 24;
    imach[6] = -125;
    imach[7] = 128;
    imach[8] = 53;
    imach[9] = -1021;
    imach[10] = 1024;

/*     MACHINE CONSTANTS FOR THE UNIVAC 1100 SERIES. */
/*
   imach[1] = 2;
   imach[2] = 35;
   imach[3] = 34359738367;
   imach[4] = 2;
   imach[5] = 27;
   imach[6] = -128;
   imach[7] = 127;
   imach[8] = 60;
   imach[9] = -1024;
   imach[10] = 1023;
*/
/*     MACHINE CONSTANTS FOR THE VAX 11/780. */
/*
   imach[1] = 2;
   imach[2] = 31;
   imach[3] = 2147483647;
   imach[4] = 2;
   imach[5] = 24;
   imach[6] = -127;
   imach[7] = 127;
   imach[8] = 56;
   imach[9] = -127;
   imach[10] = 127;
*/
    ipmpar = imach[*i];
    return ipmpar;
}


/************************************************************************
                        EVEN MORE STUFF
************************************************************************/

/*
**********************************************************************
     double genunf(double low,double high)
               GeNerate Uniform Real between LOW and HIGH
                              Function
     Generates a real uniformly distributed between LOW and HIGH.
                              Arguments
     low --> Low bound (exclusive) on real value to be generated
     high --> High bound (exclusive) on real value to be generated
**********************************************************************
*/
double genunf(double low,
              double high) {
    static double genunf;

    if (!(low > high)) goto S10;
//    cpp_printf("genunf: low > high: low=%16.6E, high=%16.6E\n", low, high);
//    _cstatfatal();
    /*NOTREACHED*/
    S10:
    genunf = low + (high - low) * ranf();
    return genunf;
}


/*
**********************************************************************
     double gengam(double a,double r)
           GENerates random deviates from GAMma distribution
                              Function
     Generates random deviates from the gamma distribution whose
     density is
          (A**R)/Gamma(R) * X**(R-1) * Exp(-A*X)
                              Arguments
     a --> Location parameter of Gamma distribution
     r --> Shape parameter of Gamma distribution
     CAREFUL: order of parameters is reversed wrt usual notation. mean is r/a; var is r/a^2
                              Method
     Renames SGAMMA from TOMS as slightly modified by BWB to use RANF
     instead of SUNIF.
     For details see:
               (Case R >= 1.0)
               Ahrens, J.H. and Dieter, U.
               Generating Gamma Variates by a
               Modified Rejection Technique.
               Comm. ACM, 25,1 (Jan. 1982), 47 - 54.
     Algorithm GD
               (Case 0.0 <= R <= 1.0)
               Ahrens, J.H. and Dieter, U.
               Computer Methods for Sampling from Gamma,
               Beta, Poisson and Binomial Distributions.
               Computing, 12 (1974), 223-246/
     Adapted algorithm GS.
**********************************************************************
*/
double gengam(double a,
              double r) {
    static double gengam;

    gengam = sgamma(r);
    gengam /= a;
    return gengam;
}


/*
**********************************************************************


     (STANDARD-)  G A M M A  DISTRIBUTION


**********************************************************************
**********************************************************************

               PARAMETER  A >= 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               GENERATING GAMMA VARIATES BY A
               MODIFIED REJECTION TECHNIQUE.
               COMM. ACM, 25,1 (JAN. 1982), 47 - 54.

     STEP NUMBERS CORRESPOND TO ALGORITHM 'GD' IN THE ABOVE PAPER
                                 (STRAIGHTFORWARD IMPLEMENTATION)

     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
     SUNIF.  The argument IR thus goes away.

**********************************************************************

               PARAMETER  0.0 < A < 1.0  !

**********************************************************************

     FOR DETAILS SEE:

               AHRENS, J.H. AND DIETER, U.
               COMPUTER METHODS FOR SAMPLING FROM GAMMA,
               BETA, POISSON AND BINOMIAL DISTRIBUTIONS.
               COMPUTING, 12 (1974), 223 - 246.

     (ADAPTED IMPLEMENTATION OF ALGORITHM 'GS' IN THE ABOVE PAPER)

**********************************************************************
     INPUT: A =PARAMETER (MEAN) OF THE STANDARD GAMMA DISTRIBUTION
     OUTPUT: SGAMMA = SAMPLE FROM THE GAMMA-(A)-DISTRIBUTION
     COEFFICIENTS Q(K) - FOR Q0 = SUM(Q(K)*A**(-K))
     COEFFICIENTS A(K) - FOR Q = Q0+(T*T/2)*SUM(A(K)*V**K)
     COEFFICIENTS E(K) - FOR EXP(Q)-1 = SUM(E(K)*Q**K)
     PREVIOUS A PRE-SET TO ZERO - AA IS A', AAA IS A"
     SQRT32 IS THE SQUAREROOT OF 32 = 5.656854249492380
*/
double sgamma(double a) {
    static double q1 = 4.166669E-2;
    static double q2 = 2.083148E-2;
    static double q3 = 8.01191E-3;
    static double q4 = 1.44121E-3;
    static double q5 = -7.388E-5;
    static double q6 = 2.4511E-4;
    static double q7 = 2.424E-4;
    static double a1 = 0.3333333;
    static double a2 = -0.250003;
    static double a3 = 0.2000062;
    static double a4 = -0.1662921;
    static double a5 = 0.1423657;
    static double a6 = -0.1367177;
    static double a7 = 0.1233795;
    static double e1 = 1.0;
    static double e2 = 0.4999897;
    static double e3 = 0.166829;
    static double e4 = 4.07753E-2;
    static double e5 = 1.0293E-2;
    static double aa = 0.0;
    static double aaa = 0.0;
    static double sqrt32 = 5.656854;
    static double sgamma, s2, s, d, t, x, u, r, q0, b, si, c, v, q, e, w, p;
    if (a == aa) goto S10;
    if (a < 1.0) goto S120;
/*
     STEP  1:  RECALCULATIONS OF S2,S,D IF A HAS CHANGED
*/
    aa = a;
    s2 = a - 0.5;
    s = sqrt(s2);
    d = sqrt32 - 12.0 * s;
    S10:
/*
     STEP  2:  T=STANDARD NORMAL DEVIATE,
               X=(S,1/2)-NORMAL DEVIATE.
               IMMEDIATE ACCEPTANCE (I)
*/
    t = snorm();
    x = s + 0.5 * t;
    sgamma = x * x;
    if (t >= 0.0) return sgamma;
/*
     STEP  3:  U= 0,1 -UNIFORM SAMPLE. SQUEEZE ACCEPTANCE (S)
*/
    u = ranf();
    if (d * u <= t * t * t) return sgamma;
/*
     STEP  4:  RECALCULATIONS OF Q0,B,SI,C IF NECESSARY
*/
    if (a == aaa) goto S40;
    aaa = a;
    r = 1.0 / a;
    q0 = ((((((q7 * r + q6) * r + q5) * r + q4) * r + q3) * r + q2) * r + q1) * r;
/*
               APPROXIMATION DEPENDING ON SIZE OF PARAMETER A
               THE CONSTANTS IN THE EXPRESSIONS FOR B, SI AND
               C WERE ESTABLISHED BY NUMERICAL EXPERIMENTS
*/
    if (a <= 3.686) goto S30;
    if (a <= 13.022) goto S20;
/*
               CASE 3:  A .GT. 13.022
*/
    b = 1.77;
    si = 0.75;
    c = 0.1515 / s;
    goto S40;
    S20:
/*
               CASE 2:  3.686 .LT. A .LE. 13.022
*/
    b = 1.654 + 7.6E-3 * s2;
    si = 1.68 / s + 0.275;
    c = 6.2E-2 / s + 2.4E-2;
    goto S40;
    S30:
/*
               CASE 1:  A .LE. 3.686
*/
    b = 0.463 + s - 0.178 * s2;
    si = 1.235;
    c = 0.195 / s - 7.9E-2 + 1.6E-2 * s;
    S40:
/*
     STEP  5:  NO QUOTIENT TEST IF X NOT POSITIVE
*/
    if (x <= 0.0) goto S70;
/*
     STEP  6:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t / (s + s);
    if (fabs(v) <= 0.25) goto S50;
    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    goto S60;
    S50:
    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
    S60:
/*
     STEP  7:  QUOTIENT ACCEPTANCE (Q)
*/
    if (log(1.0 - u) <= q) return sgamma;
    S70:
/*
     STEP  8:  E=STANDARD EXPONENTIAL DEVIATE
               U= 0,1 -UNIFORM DEVIATE
               T=(B,SI)-DOUBLE EXPONENTIAL (LAPLACE) SAMPLE
*/
    e = sexpo();
    u = ranf();
    u += (u - 1.0);
    t = b + fsign(si * e, u);
/*
     STEP  9:  REJECTION IF T .LT. TAU(1) = -.71874483771719
*/
    if (t < -0.7187449) goto S70;
/*
     STEP 10:  CALCULATION OF V AND QUOTIENT Q
*/
    v = t / (s + s);
    if (fabs(v) <= 0.25) goto S80;
    q = q0 - s * t + 0.25 * t * t + (s2 + s2) * log(1.0 + v);
    goto S90;
    S80:
    q = q0 + 0.5 * t * t * ((((((a7 * v + a6) * v + a5) * v + a4) * v + a3) * v + a2) * v + a1) * v;
    S90:
/*
     STEP 11:  HAT ACCEPTANCE (H) (IF Q NOT POSITIVE GO TO STEP 8)
*/
    if (q <= 0.0) goto S70;
    if (q <= 0.5) goto S100;
    w = exp(q) - 1.0;
    goto S110;
    S100:
    w = ((((e5 * q + e4) * q + e3) * q + e2) * q + e1) * q;
    S110:
/*
               IF T IS REJECTED, SAMPLE AGAIN AT STEP 8
*/
    if (c * fabs(u) > w * exp(e - 0.5 * t * t)) goto S70;
    x = s + 0.5 * t;
    sgamma = x * x;
    return sgamma;
    S120:
/*
     ALTERNATE METHOD FOR PARAMETERS A BELOW 1  (.3678794=EXP(-1.))
*/
    aa = 0.0;
    b = 1.0 + 0.3678794 * a;
    S130:
    p = b * ranf();
    if (p >= 1.0) goto S140;
    sgamma = exp(log(p) / a);
    if (sexpo() < sgamma) goto S130;
    return sgamma;
    S140:
    sgamma = -log((b - p) / a);
    if (sexpo() < (1.0 - a) * log(sgamma)) goto S130;
    return sgamma;
}

//
//
///*
//**********************************************************************
//
//
//     (STANDARD-)  N O R M A L  DISTRIBUTION
//
//
//**********************************************************************
//**********************************************************************
//
//     FOR DETAILS SEE:
//
//               AHRENS, J.H. AND DIETER, U.
//               EXTENSIONS OF FORSYTHE'S METHOD FOR RANDOM
//               SAMPLING FROM THE NORMAL DISTRIBUTION.
//               MATH. COMPUT., 27,124 (OCT. 1973), 927 - 937.
//
//     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM 'FL'
//     (M=5) IN THE ABOVE PAPER     (SLIGHTLY MODIFIED IMPLEMENTATION)
//
//     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
//     SUNIF.  The argument IR thus goes away.
//
//**********************************************************************
//     THE DEFINITIONS OF THE CONSTANTS A(K), D(K), T(K) AND
//     H(K) ARE ACCORDING TO THE ABOVEMENTIONED ARTICLE
//*/
double snorm(void) {
    static double a[32] = {
            0.0, 3.917609E-2, 7.841241E-2, 0.11777, 0.1573107, 0.1970991, 0.2372021, 0.2776904,
            0.3186394, 0.36013, 0.4022501, 0.4450965, 0.4887764, 0.5334097, 0.5791322,
            0.626099, 0.6744898, 0.7245144, 0.7764218, 0.8305109, 0.8871466, 0.9467818,
            1.00999, 1.077516, 1.150349, 1.229859, 1.318011, 1.417797, 1.534121, 1.67594,
            1.862732, 2.153875
    };
    static double d[31] = {
            0.0, 0.0, 0.0, 0.0, 0.0, 0.2636843, 0.2425085, 0.2255674, 0.2116342, 0.1999243,
            0.1899108, 0.1812252, 0.1736014, 0.1668419, 0.1607967, 0.1553497, 0.1504094,
            0.1459026, 0.14177, 0.1379632, 0.1344418, 0.1311722, 0.128126, 0.1252791,
            0.1226109, 0.1201036, 0.1177417, 0.1155119, 0.1134023, 0.1114027, 0.1095039
    };
    static double t[31] = {
            7.673828E-4, 2.30687E-3, 3.860618E-3, 5.438454E-3, 7.0507E-3, 8.708396E-3,
            1.042357E-2, 1.220953E-2, 1.408125E-2, 1.605579E-2, 1.81529E-2, 2.039573E-2,
            2.281177E-2, 2.543407E-2, 2.830296E-2, 3.146822E-2, 3.499233E-2, 3.895483E-2,
            4.345878E-2, 4.864035E-2, 5.468334E-2, 6.184222E-2, 7.047983E-2, 8.113195E-2,
            9.462444E-2, 0.1123001, 0.136498, 0.1716886, 0.2276241, 0.330498, 0.5847031
    };
    static double h[31] = {
            3.920617E-2, 3.932705E-2, 3.951E-2, 3.975703E-2, 4.007093E-2, 4.045533E-2,
            4.091481E-2, 4.145507E-2, 4.208311E-2, 4.280748E-2, 4.363863E-2, 4.458932E-2,
            4.567523E-2, 4.691571E-2, 4.833487E-2, 4.996298E-2, 5.183859E-2, 5.401138E-2,
            5.654656E-2, 5.95313E-2, 6.308489E-2, 6.737503E-2, 7.264544E-2, 7.926471E-2,
            8.781922E-2, 9.930398E-2, 0.11556, 0.1404344, 0.1836142, 0.2790016, 0.7010474
    };
    static long i;
    static double snorm, u, s, ustar, aa, w, y, tt;
    u = ranf();
    s = 0.0;
    if (u > 0.5) s = 1.0;
    u += (u - s);
    u = 32.0 * u;
    i = (long) (u);
    if (i == 32) i = 31;
    if (i == 0) goto S100;
/*
                                START CENTER
*/
    ustar = u - (double) i;
    aa = *(a + i - 1);
    S40:
    if (ustar <= *(t + i - 1)) goto S60;
    w = (ustar - *(t + i - 1)) * *(h + i - 1);
    S50:
/*
                                EXIT   (BOTH CASES)
*/
    y = aa + w;
    snorm = y;
    if (s == 1.0) snorm = -y;
    return snorm;
    S60:
/*
                                CENTER CONTINUED
*/
    u = ranf();
    w = u * (*(a + i) - aa);
    tt = (0.5 * w + aa) * w;
    goto S80;
    S70:
    tt = u;
    ustar = ranf();
    S80:
    if (ustar > tt) goto S50;
    u = ranf();
    if (ustar >= u) goto S70;
    ustar = ranf();
    goto S40;
    S100:
/*
                                START TAIL
*/
    i = 6;
    aa = *(a + 31);
    goto S120;
    S110:
    aa += *(d + i - 1);
    i += 1;
    S120:
    u += u;
    if (u < 1.0) goto S110;
    u -= 1.0;
    S140:
    w = u * *(d + i - 1);
    tt = (0.5 * w + aa) * w;
    goto S160;
    S150:
    tt = u;
    S160:
    ustar = ranf();
    if (ustar > tt) goto S50;
    u = ranf();
    if (ustar >= u) goto S150;
    u = ranf();
    goto S140;
}

//
//
///* Transfers sign of argument sign to argument num */
double fsign(double num,
             double sign) {
    return ((sign > 0.0 && num < 0.0) ||
            (sign < 0.0 && num > 0.0)) ? -num : num;
}

//
//
///*
//**********************************************************************
//
//
//     (STANDARD-)  E X P O N E N T I A L   DISTRIBUTION
//
//
//**********************************************************************
//**********************************************************************
//
//     FOR DETAILS SEE:
//
//               AHRENS, J.H. AND DIETER, U.
//               COMPUTER METHODS FOR SAMPLING FROM THE
//               EXPONENTIAL AND NORMAL DISTRIBUTIONS.
//               COMM. ACM, 15,10 (OCT. 1972), 873 - 882.
//
//     ALL STATEMENT NUMBERS CORRESPOND TO THE STEPS OF ALGORITHM
//     'SA' IN THE ABOVE PAPER (SLIGHTLY MODIFIED IMPLEMENTATION)
//
//     Modified by Barry W. Brown, Feb 3, 1988 to use RANF instead of
//     SUNIF.  The argument IR thus goes away.
//
//**********************************************************************
//     Q(N) = SUM(ALOG(2.0)**K/K!)    K=1,..,N ,      THE HIGHEST N
//     (HERE 8) IS DETERMINED BY Q(N)=1.0 WITHIN STANDARD PRECISION
//*/
double sexpo(void) {
    static double q[8] = {
            0.6931472, 0.9333737, 0.9888778, 0.9984959, 0.9998293, 0.9999833, 0.9999986, 1.0
    };
    static long i;
    static double sexpo, a, u, ustar, umin;
    static double *q1 = q;
    a = 0.0;
    u = ranf();
    goto S30;
    S20:
    a += *q1;
    S30:
    u += u;
    if (u <= 1.0) goto S20;
    u -= 1.0;
    if (u > *q1) goto S60;
    sexpo = a + u;
    return sexpo;
    S60:
    i = 1;
    ustar = ranf();
    umin = ustar;
    S70:
    ustar = ranf();
    if (ustar < umin) umin = ustar;
    i += 1;
    if (u > *(q + i - 1)) goto S70;
    sexpo = a + umin * *q1;
    return sexpo;
}


/*
**********************************************************************
     long mltmod(long a,long s,long m)
                    Returns (A*S) MOD M
     This is a transcription from Pascal to Fortran of routine
     MULtMod_Decompos from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     a, s, m  -->
**********************************************************************
*/
long mltmod(long a,
            long s,
            long m) {
#define h 32768L
    static long mltmod, a0, a1, k, p, q, qh, rh;
/*
     H = 2**((b-2)/2) where b = 32 because we are using a 32 bit
      machine. On a different machine recompute H
*/
    if (!(a <= 0 || a >= m || s <= 0 || s >= m)) goto S10;
//    cpp_printf("mltmod: requires (0 < a < m); (0 < s < m): ");
//    cpp_printf("a = %12ld, s = %12ld, m = %12ld\n", a, s, m);
//    _cstatfatal();
    /*NOTREACHED*/
    S10:
    if (!(a < h)) goto S20;
    a0 = a;
    p = 0;
    goto S120;
    S20:
    a1 = a / h;
    a0 = a - h * a1;
    qh = m / h;
    rh = m - h * qh;
    if (!(a1 >= h)) goto S50;
    a1 -= h;
    k = s / qh;
    p = h * (s - k * qh) - k * rh;
    S30:
    if (!(p < 0)) goto S40;
    p += m;
    goto S30;
    S40:
    goto S60;
    S50:
    p = 0;
    S60:
/*
     P = (A2*S*H)MOD M
*/
    if (!(a1 != 0)) goto S90;
    q = m / a1;
    k = s / q;
    p -= (k * (m - a1 * q));
    if (p > 0) p -= m;
    p += (a1 * (s - k * q));
    S70:
    if (!(p < 0)) goto S80;
    p += m;
    goto S70;
    S90:
    S80:
    k = p / qh;
/*
     P = ((A2*H + A1)*S)MOD M
*/
    p = h * (p - k * qh) - k * rh;
    S100:
    if (!(p < 0)) goto S110;
    p += m;
    goto S100;
    S120:
    S110:
    if (!(a0 != 0)) goto S150;
/*
     P = ((A2*H + A1)*H*S)MOD M
*/
    q = m / a0;
    k = s / q;
    p -= (k * (m - a0 * q));
    if (p > 0) p -= m;
    p += (a0 * (s - k * q));
    S130:
    if (!(p < 0)) goto S140;
    p += m;
    goto S130;
    S150:
    S140:
    mltmod = p;
    return mltmod;
#undef h
}


/*
**********************************************************************
     double ranf(void)
                RANDom number generator as a Function
     Returns a random doubleing point number from a uniform distribution
     over 0 - 1 (endpoints of this interval are not returned) using the
     current generator
     This is a transcription from Pascal to Fortran of routine
     Uniform_01 from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
double ranf(void) {
    static double ranf;
/*
     4.656613057E-10 is 1/M1  M1 is set in a data statement in IGNLGI
      and is currently 2147483563. If M1 changes, change this also.
*/
    ranf = ignlgi() * 4.656613057E-10;
    return ranf;
}


/*
**********************************************************************
     void gscgn(long getset,long *g)
                         Get/Set GeNerator
     Gets or returns in G the number of the current generator
                              Arguments
     getset --> 0 Get
                1 Set
     g <-- Number of the current random number generator (1..32)
**********************************************************************
*/
void gscgn(long getset,
           long *g) {
#define numg 32L
    static long curntg = 1;
    if (getset == 0) *g = curntg;
    else {
        if (*g < 0 || *g > numg) {
//            cpp_printf("gscgn: generator number out of range\n");
//            _cstatfatal();
            /*NOTREACHED*/
        }
        curntg = *g;
    }
#undef numg
}


/*
**********************************************************************
     void gsrgs(long getset,long *qvalue)
               Get/Set Random Generators Set
     Gets or sets whether random generators set (initialized).
     Initially (data statement) state is not set
     If getset is 1 state is set to qvalue
     If getset is 0 state returned in qvalue
**********************************************************************
*/
void gsrgs(long getset,
           long *qvalue) {
    static long qinit = 0;

    if (getset == 0) *qvalue = qinit;
    else qinit = *qvalue;
}


/*
**********************************************************************
     void gssst(long getset,long *qset)
          Get or Set whether Seed is Set
     Initialize to Seed not Set
     If getset is 1 sets state to Seed Set
     If getset is 0 returns T in qset if Seed Set
     Else returns F in qset
**********************************************************************
*/
void gssst(long getset,
           long *qset) {
    static long qstate = 0;
    if (getset != 0) qstate = 1;
    else *qset = qstate;
}


/*
**********************************************************************
     void setall(long iseed1,long iseed2)
               SET ALL random number generators
     Sets the initial seed of generator 1 to ISEED1 and ISEED2. The
     initial seeds of the other generators are set accordingly, and
     all generators states are set to these seeds.
     This is a transcription from Pascal to Fortran of routine
     Set_Initial_Seed from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     iseed1 -> First of two integer seeds
     iseed2 -> Second of two integer seeds
**********************************************************************
*/
void setall(long iseed1,
            long iseed2) {
#define numg 32L
    extern void gsrgs(long getset, long *qvalue);
    extern void gssst(long getset, long *qset);
    extern void gscgn(long getset, long *g);
    extern long Xm1, Xm2, Xa1vw, Xa2vw, Xig1[32], Xig2[32];
    static long T1;
    static long g, ocgn;
    static long qrgnin;
    T1 = 1;
/*
     TELL IGNLGI, THE ACTUAL NUMBER GENERATOR, THAT THIS ROUTINE
      HAS BEEN CALLED.
*/
    gssst(1, &T1);
    gscgn(0L, &ocgn);
/*
     Initialize Common Block if Necessary
*/
    gsrgs(0L, &qrgnin);
    if (!qrgnin) inrgcm();
    *Xig1 = iseed1;
    *Xig2 = iseed2;
    initgn(-1L);
    for (g = 2; g <= numg; g++) {
        *(Xig1 + g - 1) = mltmod(Xa1vw, *(Xig1 + g - 2), Xm1);
        *(Xig2 + g - 1) = mltmod(Xa2vw, *(Xig2 + g - 2), Xm2);
        gscgn(1L, &g);
        initgn(-1L);
    }
    gscgn(1L, &ocgn);
#undef numg
}


/*
**********************************************************************
     void initgn(long isdtyp)
          INIT-ialize current G-e-N-erator
     Reinitializes the state of the current generator
     This is a transcription from Pascal to Fortran of routine
     Init_Generator from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
                              Arguments
     isdtyp -> The state to which the generator is to be set
          isdtyp = -1  => sets the seeds to their initial value
          isdtyp =  0  => sets the seeds to the first value of
                          the current block
          isdtyp =  1  => sets the seeds to the first value of
                          the next block
**********************************************************************
*/
void initgn(long isdtyp) {
#define numg 32L
    extern void gsrgs(long getset, long *qvalue);
    extern void gscgn(long getset, long *g);
    extern long Xm1, Xm2, Xa1w, Xa2w, Xig1[32], Xig2[32], Xlg1[32], Xlg2[32], Xcg1[32], Xcg2[32];
    static long g;
    static long qrgnin;
/*
     Abort unless random number generator initialized
*/
    gsrgs(0L, &qrgnin);
    if (qrgnin) goto S10;
//    cpp_printf("initgn: random number generator not initialized\n");
//    _cstatfatal();
    /*NOTREACHED*/
    S10:
    gscgn(0L, &g);
    if (-1 != isdtyp) goto S20;
    *(Xlg1 + g - 1) = *(Xig1 + g - 1);
    *(Xlg2 + g - 1) = *(Xig2 + g - 1);
    goto S50;
    S20:
    if (0 != isdtyp) goto S30;
    goto S50;
    S30:
/*
     do nothing
*/
    if (1 != isdtyp) goto S40;
    *(Xlg1 + g - 1) = mltmod(Xa1w, *(Xlg1 + g - 1), Xm1);
    *(Xlg2 + g - 1) = mltmod(Xa2w, *(Xlg2 + g - 1), Xm2);
    goto S50;
    S40:
//    cpp_printf("initgn: isdtyp not in range\n");
//    _cstatfatal();
    /*NOTREACHED*/
    S50:
    *(Xcg1 + g - 1) = *(Xlg1 + g - 1);
    *(Xcg2 + g - 1) = *(Xlg2 + g - 1);
#undef numg
}


/*
**********************************************************************
     long ignlgi(void)
               GeNerate LarGe Integer
     Returns a random integer following a uniform distribution over
     (1, 2147483562) using the current generator.
     This is a transcription from Pascal to Fortran of routine
     Random from the paper
     L'Ecuyer, P. and Cote, S. "Implementing a Random Number Package
     with Splitting Facilities." ACM Transactions on Mathematical
     Software, 17:98-111 (1991)
**********************************************************************
*/
long ignlgi(void) {
#define numg 32L
    extern void gsrgs(long getset, long *qvalue);
    extern void gssst(long getset, long *qset);
    extern void gscgn(long getset, long *g);
    extern long Xm1, Xm2, Xa1, Xa2, Xcg1[32], Xcg2[32];
    extern long Xqanti[32];
    static long ignlgi, curntg, k, s1, s2, z;
    static long qqssd, qrgnin;
/*
     IF THE RANDOM NUMBER PACKAGE HAS NOT BEEN INITIALIZED YET, DO SO.
     IT CAN BE INITIALIZED IN ONE OF TWO WAYS : 1) THE FIRST CALL TO
     THIS ROUTINE  2) A CALL TO SETALL.
*/
    gsrgs(0L, &qrgnin);
    if (!qrgnin) inrgcm();
    gssst(0, &qqssd);
    if (!qqssd) setall(1234567890L, 123456789L);
/*
     Get Current Generator
*/
    gscgn(0L, &curntg);
    s1 = *(Xcg1 + curntg - 1);
    s2 = *(Xcg2 + curntg - 1);
    k = s1 / 53668L;
    s1 = Xa1 * (s1 - k * 53668L) - k * 12211;
    if (s1 < 0) s1 += Xm1;
    k = s2 / 52774L;
    s2 = Xa2 * (s2 - k * 52774L) - k * 3791;
    if (s2 < 0) s2 += Xm2;
    *(Xcg1 + curntg - 1) = s1;
    *(Xcg2 + curntg - 1) = s2;
    z = s1 - s2;
    if (z < 1) z += (Xm1 - 1);
    if (*(Xqanti + curntg - 1)) z = Xm1 - z;
    ignlgi = z;
    return ignlgi;
#undef numg
}


/*
**********************************************************************
     void inrgcm(void)
          INitialize Random number Generator CoMmon
                              Function
     Initializes common area  for random number  generator.  This saves
     the  nuisance  of  a  BLOCK DATA  routine  and the  difficulty  of
     assuring that the routine is loaded with the other routines.
**********************************************************************
*/
void inrgcm(void) {
#define numg 32L
    extern void gsrgs(long getset, long *qvalue);
    extern long Xm1, Xm2, Xa1, Xa2, Xa1w, Xa2w, Xa1vw, Xa2vw;
    extern long Xqanti[32];
    static long T1;
    static long i;
/*
     V=20;                            W=30;
     A1W = MOD(A1**(2**W),M1)         A2W = MOD(A2**(2**W),M2)
     A1VW = MOD(A1**(2**(V+W)),M1)    A2VW = MOD(A2**(2**(V+W)),M2)
   If V or W is changed A1W, A2W, A1VW, and A2VW need to be recomputed.
    An efficient way to precompute a**(2*j) MOD m is to start with
    a and square it j times modulo m using the function MLTMOD.
*/
    Xm1 = 2147483563L;
    Xm2 = 2147483399L;
    Xa1 = 40014L;
    Xa2 = 40692L;
    Xa1w = 1033780774L;
    Xa2w = 1494757890L;
    Xa1vw = 2082007225L;
    Xa2vw = 784306273L;
    for (i = 0; i < numg; i++) *(Xqanti + i) = 0;
    T1 = 1;
/*
     Tell the world that common has been initialized
*/
    gsrgs(1L, &T1);
#undef numg
}


///************************************************************************
//            FUNCTION OPTIMIZATION
//************************************************************************/
//
///*
// * PURPOSE: UNIVARIATE OPTIMIZATION WITHOUT DERIVATIVE INFORMATION
// *
// * Given a function f and a bracketing triplet of abscissas ax, bx, cx
// * (such that bx is between ax & cx and f(bx) is less than f(ax) and f(cx)),
// * this routine isolates the minimum to a fractional precision of about eps
// * using Brent's method. The abscissa of the minimum is returned as xmin and
// * the value of the function at the minimum is the returned value.
// */
double univmin(double ax,
               double bx,
               double cx,
               double (*const f)(double),
               double eps,
               double *xmin,
               int itmax) {
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d);

    const double CGOLD = 0.3819660;  //golden ratio
    const double ZEPS = 1.0e-10;   //protects against trying to achieve fractional accuracy when min is exactly zero
    int iter;
    double a, b, d = 1, etemp, fu, fv, fw, fx, p, q, r, eps1, eps2, u, v, w, x, xm;
    double e = 0.0;

    //assert(f != NULL);
    //assert(xmin != NULL);

    a = (ax < cx ? ax : cx);     //a,b must be in ascending order but input abscissas need not be
    b = (ax > cx ? ax : cx);

    x = w = v = bx;                  //initializations
    fw = fv = fx = (*f)(x);

    for (iter = 1; iter <= itmax; iter++) {      //main program loop
        xm = 0.5 * (a + b);
        eps2 = 2.0 * (eps1 = eps * fabs(x) + ZEPS);
        if (fabs(x - xm) <= (eps2 - 0.5 * (b - a))) {  //test for done here
            *xmin = x;
            return fx;
        }
        if (fabs(e) > eps1) {                  //construct a trial parabolic fit
            r = (x - w) * (fx - fv);
            q = (x - v) * (fx - fw);
            p = (x - v) * q - (x - w) * r;
            q = 2.0 * (q - r);
            if (q > 0.0) p = -p;
            q = fabs(q);
            etemp = e;
            e = d;
            //determine acceptability of parabolic fit & take the golden section step into the larger of the two segments
            if (fabs(p) >= fabs(0.5 * q * etemp) || p <= q * (a - x) || p >= q * (b - x)) {
                d = CGOLD * (e = (x >= xm ? a - x : b - x));
            } else {
                d = p / q;                            //take the parabolic step
                u = x + d;
                if (u - a < eps2 || b - u < eps2) d = SETSIGN(eps1, xm - x);
            }
        } else {
            d = CGOLD * (e = (x >= xm ? a - x : b - x));
        }
        u = (fabs(d) >= eps1 ? x + d : x + SETSIGN(eps1, d));
        fu = (*f)(u);                         //this is the one function evaluation per iteration
        if (fu <= fx) {                      //now decide what to do with our function evaluation
            if (u >= x) a = x; else b = x;
            SHFT(v, w, x, u);                     //housekeeping
            SHFT(fv, fw, fx, fu);
        } else {
            if (u < x) a = u; else b = u;
            if (fu <= fw || w == x) {
                v = w;
                w = u;
                fv = fw;
                fw = fu;
            } else if (fu <= fv || v == x || v == w) {
                v = u;
                fv = fu;
            }
        }                                   //done with housekeeping. Back for another iteration
    }

    *xmin = x;                             //only get here if iteration limit is reached
    return fx;
#undef SHFT
}

void minimize(double th[],
              double **dirini,
              int n,
              double eps,
              int *iter,
              double *fret,
              double (*const f)(double []),
              int itmax) {
    int i;
    int j;
    int ibig;
    double del;
    double fth;
    double fthtt;
    double t;
    double *tht;
    double *thtt;
    double *dirinit;
    bool converged = false;

    //assert(dirini != NULL);
    //assert(iter != NULL);
    //assert(fret != NULL);
    //assert(f != NULL);

    tht = dvector(1, n);
    thtt = dvector(1, n);
    dirinit = dvector(1, n);

    /* Initial parameter and function values */
    *fret = (*f)(th);
    for (j = 1; j <= n; j++) {
        tht[j] = th[j];
    }

    for (*iter = 1; (*iter < itmax) && (converged == false); ++(*iter)) {
        fth = (*fret);
        ibig = 0;
        del = 0.0;

        /* Minimize along all directions */
        for (i = 1; i <= n; i++) {
            for (j = 1; j <= n; j++) {
                dirinit[j] = dirini[j][i];
            }
            fthtt = (*fret);
            dirmin(th, dirinit, n, fret, f, itmax, eps);
            if (fabs(fthtt - (*fret)) > del) {
                del = fabs(fthtt - (*fret));
                ibig = i;
            }
        }

        /*
         * extrapolated point, average direction and function value
         * at extrapolated point
         */
        for (j = 1; j <= n; j++) {
            thtt[j] = 2.0 * th[j] - tht[j];
            dirinit[j] = th[j] - tht[j];
            tht[j] = th[j];
        }
        fthtt = (*f)(thtt);
        if (fthtt < fth) {
            t = 2.0 * (fth - 2.0 * (*fret) + fthtt) *
                sqrt(fth - (*fret) - del) - del * sqrt(fth - fthtt);
            if (t < 0.0) {
                dirmin(th, dirinit, n, fret, f, itmax, eps);
                for (j = 1; j <= n; j++) {
                    dirini[j][ibig] = dirini[j][n];
                    dirini[j][n] = dirinit[j];
                }
            }
        }

        if (2.0 * fabs(fth - (*fret)) <= eps * (fabs(fth) + fabs(*fret))) {
            converged = true;
        }
    }

    free_dvector(dirinit, 1, n);
    free_dvector(thtt, 1, n);
    free_dvector(tht, 1, n);
}

//
//
///*
// * Global variables communicate with f1dim.
// */
// TODO OMG PLS NO GLOBAL VARS
int ncom;
double *pcom;
double *xicom;

//
double (*nrfunc)(double []);

//
//
///* Must accompany dirmin. */
double f1dim(double x) {
    int j;
    double f;
    double *xt;

    xt = dvector(1, ncom);
    for (j = 1; j <= ncom; j++) {
        xt[j] = pcom[j] + x * xicom[j];
    }
    f = (*nrfunc)(xt);
    free_dvector(xt, 1, ncom);

    return f;
}

//
//
///*
// * Given an n-dimensional point p[1..n] and an n-dimensional direction xi[1..n],
// * moves and resets p to where the function func(p) takes on a minimum along
// * the direction xi from p, and replaces xi by the actual vector displacement
// * that p was moved. Also returns as fret the value of func at the returned
// * location p. This is actually all accomplished by calling the routines
// * mnbrak() and univmin().
// */
void dirmin(double p[],
            double xi[],
            int n,
            double *fret,
            double (*const func)(double []),
            int itmax,
            double dirminEPS) {
    int j;
    double xmin;
    double fx;
    double fb;
    double fa;
    double bx;
    double ax = 0.0;
    double xx = 1.0;

    //assert(fret != NULL);
    //assert(func != NULL);

    /* Set the global variables */
    ncom = n;
    pcom = dvector(1, n);
    xicom = dvector(1, n);
    nrfunc = func;
    for (j = 1; j <= n; j++) {
        pcom[j] = p[j];
        xicom[j] = xi[j];
    }

    mnbrak(&ax, &xx, &bx, &fa, &fx, &fb, f1dim);
    *fret = univmin(ax, xx, bx, f1dim, dirminEPS, &xmin, itmax);

    /* Construct the vector results to return */
    for (j = 1; j <= n; j++) {
        xi[j] *= xmin;
        p[j] += xi[j];
    }
    free_dvector(xicom, 1, n);
    free_dvector(pcom, 1, n);
}

//
//
static double maxarg1;
static double maxarg2;

#define FMAX(a, b) (maxarg1=(a),maxarg2=(b),(maxarg1) > (maxarg2) ? (maxarg1) : (maxarg2))

//
///*
// * Given a function func, and given distinct initial points ax and bx, search
// * in the downhill direction (defined by the function as evaluated at the
// * initial points) and returns new points ax, bx, cx that bracket a minimum
// * of the function. Also returned are the function values at the three points,
// * fa, fb, and fc.
// */
void mnbrak(double *ax,
            double *bx,
            double *cx,
            double *fa,
            double *fb,
            double *fc,
            double (*const func)(double)) {
#define SHFT(a, b, c, d) (a)=(b);(b)=(c);(c)=(d)

    const double GOLD = 1.618034;  /* default ratio by which successive intervals are magnified */
    const double GLIMIT = 100.0;   /* maximum magnification allowed for a parabolic-fit step */
    const double TINY = 1.0e-25;

    //assert(ax != NULL);
    //assert(bx != NULL);
    //assert(cx != NULL);
    //assert(fa != NULL);
    //assert(fb != NULL);
    //assert(fc != NULL);
    //assert(func != NULL);

    *fa = (*func)(*ax);
    *fb = (*func)(*bx);
    if (*fb > *fa) {
        double dum;

        /*
         * Switch roles of a and b so that we can go downhill
         * in the direction from a to b.
         */
        SHFT(dum, *ax, *bx, dum);
        SHFT(dum, *fb, *fa, dum);
    }
    *cx = (*bx) + GOLD * (*bx - *ax);  /* First guess for c */
    *fc = (*func)(*cx);

    /* Keep returning here until we bracket... */
    while (*fb > *fc) {
        double r;
        double q;
        double u;
        double ulim;
        double fu;

        /*
         * Compute u by parabolic extrapolation from a, b, c.
         * TINY prevents any division by zero.
         */
        r = (*bx - *ax) * (*fb - *fc);
        q = (*bx - *cx) * (*fb - *fa);
        u = (*bx) - ((*bx - *cx) * q - (*bx - *ax) * r) /
                    (2.0 * SETSIGN(FMAX(fabs(q - r), TINY), q - r));
        ulim = (*bx) + GLIMIT * (*cx - *bx); /* Go no farther than this. */

        /* Test various possibilities */
        if ((*bx - u) * (u - *cx) > 0.0) {
            /* Parabolic u is between b and c: try it */
            fu = (*func)(u);
            if (fu < *fc) {
                /* Got a minimum between b and c */
                *ax = (*bx);
                *bx = u;
                *fa = (*fb);
                *fb = fu;
                return;
            } else if (fu > *fb) {
                /* Got a minimum between between a and u */
                *cx = u;
                *fc = fu;
                return;
            }
            /* Parabolic fit was no use, use default magnification */
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        } else if ((*cx - u) * (u - ulim) > 0.0) {
            /* Parabolic fit is between c and its allowed limit */
            fu = (*func)(u);
            if (fu < *fc) {
                SHFT(*bx, *cx, u, *cx + GOLD * (*cx - *bx));
                SHFT(*fb, *fc, fu, (*func)(u));
            }
        } else if ((u - ulim) * (ulim - *cx) >= 0.0) {
            /* Limit parabolic u to maximum allowed value */
            u = ulim;
            fu = (*func)(u);
        } else {
            /* Reject parabolic u, use default magnification */
            u = (*cx) + GOLD * (*cx - *bx);
            fu = (*func)(u);
        }

        /* Eliminate oldest point and continue */
        SHFT(*ax, *bx, *cx, u);
        SHFT(*fa, *fb, *fc, fu);
    }
#undef SHFT
}
//
