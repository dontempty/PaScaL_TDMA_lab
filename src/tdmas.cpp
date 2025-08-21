#include <vector>
#include <iostream>
#include "tdmas.hpp"

void tdma_single(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n1) {

    // @brief       Solve a single tridiagonal system of equations using the Thomas algorithm.
    // @param       a       Coefficients in lower diagonal elements
    // @param       b       Coefficients in diagonal elements
    // @param       c       Coefficients in upper diagonal elements
    // @param       d       Coefficients in the right-hand side terms
    // @param       n1      Number of rows in each process, dimension of tridiagonal matrix N divided by nprocs

    int i;
    double r;

    d[0] = d[0]/b[0];
    c[0] = c[0]/b[0];

    for (i=1; i<n1; ++i) {
        r = 1.0/(b[i]-a[i]*c[i-1]);
        d[i] = r*(d[i]-a[i]*d[i-1]);
        c[i] = r*c[i];
    }

    for (i=n1-2; i>=0; --i) {
        d[i] = d[i]-c[i]*d[i+1];
    }
}

void tdma_cycl_single(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n1) {

    // @brief       Solve a single cyclic tridiagonal system of equations using the Thomas algorithm.
    // @param       a       Coefficients in lower diagonal elements
    // @param       b       Coefficients in diagonal elements
    // @param       c       Coefficients in upper diagonal elements
    // @param       d       Coefficients in the right-hand side terms
    // @param       n1      Number of rows in each process, dimension of tridiagonal matrix N divided by nprocs
    
    int i;
    double rr;

    std::vector<double> e(n1, 0.0);
    e[1] = -a[1];
    e[n1-1] = -c[n1-1];

    d[1] = d[1]/b[1];
    e[1] = e[1]/b[1];
    c[1] = c[1]/b[1];

    for (i=2; i<=n1-1; ++i) {
        rr = 1.0/(b[i]-a[i]*c[i-1]);
        d[i] = rr*(d[i]-a[i]*d[i-1]);
        e[i] = rr*(e[i]-a[i]*e[i-1]);
        c[i] = rr*c[i];
    }

    for (i=n1-2; i>=1; --i) {
        d[i] = d[i]-c[i]*d[i+1];
        e[i] = e[i]-c[i]*e[i+1];
    }

    d[0] = (d[0] - a[0]*d[n1-1] - c[0]*d[1])/(b[0] + a[0]*e[n1-1] + c[0]*e[1]);

    for (i=1; i<=n1-1; ++i) {
        d[i] = d[i] + d[0]*e[i];
    }
}

void tdma_many(
    std::vector<double> &a,
    std::vector<double> &b,
    std::vector<double> &c,
    std::vector<double> &d,
    int n1, int n2) {
    
    // @brief       Solve many tridiagonal systems of equations using the Thomas algorithm.
    //              First index indicates the number of independent many tridiagonal systems to use vectorization.
    //              Second index indicates the row number in the tridiagonal system .
    // @param       a       Coefficient array in lower diagonal elements
    // @param       b       Coefficient array in diagonal elements
    // @param       c       Coefficient array in upper diagonal elements
    // @param       d       Coefficient array in the right-hand side terms
    // @param       n1      Number of tridiagonal systems per process
    // @param       n2      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs
        
    int idx;
    int i, j;
    double r;

    // Forward elimination
    for (j=0; j<n1; ++j) {
        idx = j*n2 + 0;

        d[idx] /= b[idx];
        c[idx] /= b[idx];
    }

    for (j=0; j<n1; ++j) {
        for (i = 1; i < n2; ++i) {
            idx = j*n2 + i;
            r = 1.0 / (b[idx] - a[idx] * c[idx - 1]);
            d[idx] = r * (d[idx] - a[idx] * d[idx - 1]);
            c[idx] = r * c[idx];
        }
    }

    // Back substitution
    for (j=0; j<n1; ++j) {
        for (i=n2-2; i>=0; --i) {
            idx = j*n2 + i;

            d[idx] = d[idx] - c[idx] * d[idx + 1];
        }
    }
}

void tdma_cycl_many(
    std::vector<double> &a,
    std::vector<double> &b,
    std::vector<double> &c,
    std::vector<double> &d,
    int n1, int n2) {
    
    // @brief       Solve many cyclic tridiagonal systems of equations using the Thomas algorithm.
    //              First index indicates the number of independent many tridiagonal systems to use vectorization.
    //              Second index indicates the row number in the tridiagonal system.
    // @param       a       Coefficient array in lower diagonal elements
    // @param       b       Coefficient array in diagonal elements
    // @param       c       Coefficient array in upper diagonal elements
    // @param       d       Coefficient array in the right-hand side terms
    // @param       n1      Number of tridiagonal systems per process
    // @param       n2      Number of rows in each process, size of the tridiagonal matrix N divided by nprocs

    std::vector<double> e(n1*n2);
    int idx;
    int i, j;
    double r;

    for (j=0; j<n1; ++j) {
        for (i=0; i<n2; ++i) {
            idx = j*n2 +i;
            e[idx] = 0;
        }
        idx = j*n2 + 1;
        e[idx] = -a[idx];
        
        idx = j*n2 + (n2-1);
        e[idx] = -c[idx];
    }

    for (j=0; j<n1; ++j) {
        idx = j*n2 + 1;
        d[idx] /= b[idx];
        e[idx] /= b[idx];
        c[idx] /= b[idx];
    }

    for (j=0; j<n1; ++j) {
        for (i=2; i<n2; ++i) {
            idx = j*n2 + i;
            r = 1.0 / (b[idx] - a[idx]*c[idx-1]);
            d[idx] = r * (d[idx] - a[idx]*d[idx-1]);
            e[idx] = r * (e[idx] - a[idx]*e[idx-1]);
            c[idx] = r * c[idx];
        }
    }

    for (j=0; j<n1; ++j) {
        for (i=n2-2; i>=1; --i) {
            idx = j*n2 + i;
            d[idx] = d[idx] - c[idx] * d[idx + 1];
            e[idx] = e[idx] - c[idx] * e[idx + 1];
        }
    }

    for (j=0; j<n1; ++j) {
        idx = j*n2 + 0;
        d[idx] = (d[idx] - a[idx]*d[idx + (n2-1)] - c[idx]*d[idx + 1]) \
               / (b[idx] + a[idx]*e[idx + (n2-1)] + c[idx]*e[idx + 1]);
    }

    double dd;
    for (j=0; j<n1; ++j) {
        dd = d[j*n2 + 0];
        for (i=1; i<n2; ++i) {
            idx = j*n2 + i;
            d[idx] = d[idx] + dd*e[idx];
        }
    }
}

