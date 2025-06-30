#include <vector>
#include <iostream>
#include "tdmas.hpp"

void tdma_single(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n1) {

    int i;
    double r;

    d[0] = d[0]/b[0];
    c[0] = c[0]/b[0];

    for (i=1; i<=n1-1; ++i) {
        r = 1.0/(b[i]-a[i]*c[i-1]);
        d[i] = r*(d[i]-a[i]*d[i-1]);
        c[i] = r*c[i];
    }

    for (i=n1-2; i>=0; --i) {
        d[i] = d[i]-c[i]*d[i+1];
    }
}

void tdma_cycl_single(std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d, int n1) {

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

    d[0] = (d[0]-a[0]*d[n1-1]-c[0]*d[1])/(b[0]+a[0]*e[n1-1]+c[0]*e[1]);

    for (i=1; i<=n1-1; ++i) {
        d[i] = d[i] + d[0]*e[i];
    }
}

void tdma_many(
    std::vector<std::vector<double>> &a,
    std::vector<std::vector<double>> &b,
    std::vector<std::vector<double>> &c,
    std::vector<std::vector<double>> &d,
    int n1, int n2) {
        
    std::vector<double> r(n1);

    // Forward elimination
    for (int i = 0; i < n1; ++i) {
        d[i][0] /= b[i][0];
        c[i][0] /= b[i][0];
    }

    for (int j = 1; j < n2; ++j) {
        for (int i = 0; i < n1; ++i) {
            r[i] = 1.0 / (b[i][j] - a[i][j] * c[i][j - 1]);
            d[i][j] = r[i] * (d[i][j] - a[i][j] * d[i][j - 1]);
            c[i][j] = r[i] * c[i][j];
        }
    }

    // Back substitution
    for (int j = n2 - 2; j >= 0; --j) {
        for (int i = 0; i < n1; ++i) {
            d[i][j] = d[i][j] - c[i][j] * d[i][j + 1];
        }
    }
}

