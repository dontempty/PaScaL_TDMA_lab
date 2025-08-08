#include "../examples_lab/save.hpp"
#include "../examples_lab/index.hpp"

#include "iostream"
#include <vector>
#include <cmath>
#include <chrono> 

void tdma_many(
    std::vector<double> &a,
    std::vector<double> &b,
    std::vector<double> &c,
    std::vector<double> &d,
    int n1, int n2) {
    // n1: n_sys
    // n2: n_row
        
    std::vector<double> r(n1);
    int idx;
    int i, j;

    // Forward elimination
    for (j = 0; j < n1; ++j) {
        idx = j*n2 + 0;

        d[idx] /= b[idx];
        c[idx] /= b[idx];
    }

    for (j = 0; j < n1; ++j) {
        // r[j] = 1.0 / (b[idx] - a[idx] * c[idx - 1]);
        for (i = 1; i < n2; ++i) {
            idx = j*n2 + i;
            double r = 1.0 / (b[idx] - a[idx] * c[idx - 1]);
            d[idx] = r * (d[idx] - a[idx] * d[idx - 1]);
            c[idx] = r * c[idx];
        }
    }

    // Back substitution
    for (j = 0; j < n1; ++j) {
        for (i = n2-2; i >= 0; --i) {
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
    // n1: n_sys
    // n2: n_row
        
    // std::vector<double> r(n2);
    std::vector<double> e(n1*n2);
    int idx;
    int i, j;

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

    double r;
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

// compile = g++ -O2 -std=c++17 -static-libgcc -static-libstdc++ -o many_fixed_P_3D.exe many_fixed_P_3D.cpp
// Run = ./many_fixed_P_3D.exe 
int main() {

    int i, j, k;

    // paramater
    int Nx = 128;
    int Ny = 128;
    int Nz = 128;
    int nx1 = Nx+2;
    int ny1 = Ny+2;
    int nz1 = Nz+2;

    double x0 = -1;
    double xN = 1;
    double y0 = -1;
    double yN = 1;
    double z0 = -1;
    double zN = 1;

    double dy = (yN-y0)/Ny;
    double dx = (xN-x0)/Nx;
    double dz = (zN-z0)/Nz;

    std::vector<double> X(nx1);
    std::vector<double> Y(ny1);
    std::vector<double> Z(nz1);
    std::vector<double> DX(nx1);
    std::vector<double> DY(ny1);
    std::vector<double> DZ(nz1);
    for (i=0; i<nx1; ++i) {
        if (i==0) {
            X[i] = x0;
            Y[i] = y0;
            Z[i] = z0;
        }
        else if (i==nx1-1) {
            X[i] = xN;
            Y[i] = yN;
            Z[i] = zN;
        }
        else {
            X[i] = x0 + dx/2 + (i-1)*dx;
            Y[i] = y0 + dy/2 + (i-1)*dy;
            Z[i] = z0 + dz/2 + (i-1)*dz;
        }
        DX[i] = dx;
        DY[i] = dy;
        DZ[i] = dz;
    }

    std::vector<double> theta_x_left_index(nx1, 0.0);
    std::vector<double> theta_x_right_index(nx1, 0.0);
    std::vector<double> theta_y_left_index(ny1, 0.0);
    std::vector<double> theta_y_right_index(ny1, 0.0);
    std::vector<double> theta_z_left_index(nz1, 0.0);
    std::vector<double> theta_z_right_index(nz1, 0.0);
    theta_x_left_index[1] = 1;
    theta_x_right_index[nx1-2] = 1;
    theta_y_left_index[1] = 1;
    theta_y_right_index[ny1-2] = 1;
    theta_z_left_index[1] = 1;
    theta_z_right_index[nz1-2] = 1;

    std::vector<double> Axx((nx1-2)*(ny1-2)), Bxx((nx1-2)*(ny1-2)), Cxx((nx1-2)*(ny1-2)), Dxx((nx1-2)*(ny1-2));
    std::vector<double> Ayy((ny1-2)*(nz1-2)), Byy((ny1-2)*(nz1-2)), Cyy((ny1-2)*(nz1-2)), Dyy((ny1-2)*(nz1-2));
    std::vector<double> Azz((nz1-2)*(nx1-2)), Bzz((nz1-2)*(nx1-2)), Czz((nz1-2)*(nx1-2)), Dzz((nz1-2)*(nx1-2));

    std::vector<double> rhs_x(nx1 * ny1 * nz1, 0.0);
    std::vector<double> rhs_y(nx1 * ny1 * nz1, 0.0);
    std::vector<double> rhs_z(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta_z(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta_y(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta(nx1 * ny1 * nz1, 0.0);
    
    int idx;
    int ij, jk, ik, ki;
    int ijk, jki, kij;
    int idx_ip, idx_im;
    int idx_jp, idx_jm;
    int idx_kp, idx_km;
    double dxdx, dydy, dzdz;
    double coef_x_a, coef_x_b, coef_x_c;
    double coef_y_a, coef_y_b, coef_y_c;
    double coef_z_a, coef_z_b, coef_z_c;
    double Pi = 3.14159265358979323846;

    // theta Init
    for (k=0; k<nz1; ++k) {
        for (j=0; j<ny1; ++j) {
            for (i=0; i<nx1; ++i) {
                ijk = idx_ijk(i, j, k, nx1, ny1);

                theta[ijk] = sin(Pi*X[i]) * sin(Pi*Y[j]) * sin(Pi*Z[k]) * exp(-3.0 * Pi*Pi * 0.0 ) + cos(Pi*X[i]) * cos(Pi*Y[j]) * cos(Pi*Z[k]);
            }
        }
    }

    // bdy condition (x, y, z) = (periodic, Dirichlet, periodic)

    // bdy variable
    std::vector<double> x_left_bdy(ny1 * nz1, 0.0);
    std::vector<double> x_right_bdy(ny1 * nz1, 0.0);
    for (k=0; k<nz1; ++k) {
        for (j=0; j<ny1; ++j) {
            jk = idx_jk(j, k, ny1);

            ijk = idx_ijk(nx1-2, j, k, nx1, ny1);
            x_left_bdy[jk] = theta[ijk];

            ijk = idx_ijk(1, j, k, nx1, ny1);
            x_right_bdy[jk] = theta[ijk];
        }
    }

    std::vector<double> y_left_bdy(nz1 * nx1, 0.0);
    std::vector<double> y_right_bdy(nz1 * nx1, 0.0);
    for (k=0; k<nz1; ++k) {
        for (i=0; i<nx1; ++i) {
            ik = idx_ik(i, k, nx1);

            ijk = idx_ijk(i, 0, k, nx1, ny1);
            y_left_bdy[ik] = theta[ijk];

            ijk = idx_ijk(i, ny1-1, k, nx1, ny1);
            y_right_bdy[ik] = theta[ijk];
        }
    }

    std::vector<double> z_left_bdy(nx1 * ny1, 0.0);
    std::vector<double> z_right_bdy(nx1 * ny1, 0.0);
    for (j=0; j<ny1; ++j) {
        for (i=0; i<nx1; ++i) {
            ij = idx_ij(i, j, nx1);

            ijk = idx_ijk(i, j, nz1-2, nx1, ny1);
            z_left_bdy[ij] = theta[ijk];

            ijk = idx_ijk(i, j, 1, nx1, ny1);
            z_right_bdy[ij] = theta[ijk];
        }
    }

    int max_iter = 10;
    double dt = 0.0005;
    int time;
    auto start = std::chrono::steady_clock::now();
    for (time=0; time<max_iter; ++time) {

        // update ghost cell

        // x
        for (k=0; k<nz1; ++k) {
            for (j=0; j<ny1; ++j) {

                ijk = idx_ijk(nx1-2, j, k, nx1, ny1);
                theta[idx_ijk(0, j, k, nx1, ny1)] = theta[ijk];

                ijk = idx_ijk(1, j, k, nx1, ny1);
                theta[idx_ijk(nx1-1, j, k, nx1, ny1)] = theta[ijk];
            }
        }

        // y
        for (k=0; k<nz1; ++k) {
            for (i=0; i<nx1; ++i) {
                ik = idx_ik(i, k, nx1);

                ijk = idx_ijk(i, 0, k, nx1, ny1);
                y_left_bdy[ik] = theta[ijk];

                ijk = idx_ijk(i, ny1-1, k, nx1, ny1);
                y_right_bdy[ik] = theta[ijk];
            }
        }


        // z
        for (j=0; j<ny1; ++j) {
            for (i=0; i<nx1; ++i) {

                ijk = idx_ijk(i, j, nz1-2, nx1, ny1);
                theta[idx_ijk(i, j, 0, nx1, ny1)] = theta[ijk];

                ijk = idx_ijk(i, j, 1, nx1, ny1);
                theta[idx_ijk(i, j, nz1-1, nx1, ny1)] = theta[ijk];
            }
        }

        // Calculating r.h.s -------------------------------------------------------------------

        // rhs_x ---------------------------
        for (k=0; k<nz1; ++k) {
            for (j=0; j<ny1; ++j) {
                for (i=1; i<nx1-1; ++i) {
                    ijk    = idx_ijk(i  , j, k, nx1, ny1);
                    idx_ip = idx_ijk(i+1, j, k, nx1, ny1);
                    idx_im = idx_ijk(i-1, j, k, nx1, ny1);
                    dxdx = DX[i]*DX[i];

                    coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 );
                    coef_x_b = (dt / 2.0 / dxdx) * (-2.0 );
                    coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 );

                    rhs_x[ijk] = (coef_x_c*theta[idx_ip] + (1.0+coef_x_b)*theta[ijk] + coef_x_a*theta[idx_im]);
                }
            }
        }

        // rhs_y ---------------------------
        for (i=1; i<nx1-1; ++i) {
            for (k=0; k<nz1; ++k) {
                for (j=1; j<ny1-1; ++j) {
                    ijk    = idx_ijk(i, j  , k, nx1, ny1);
                    idx_jp = idx_ijk(i, j+1, k, nx1, ny1);
                    idx_jm = idx_ijk(i, j-1, k, nx1, ny1);
                    dydy = DY[j]*DY[j];

                    coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * theta_y_left_index[j] + (1.0/3.0) * theta_y_right_index[j] );
                    coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -     (2.0) * theta_y_right_index[j] );
                    coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] + (5.0/3.0) * theta_y_right_index[j] );

                    rhs_y[ijk] = (coef_y_c*rhs_x[idx_jp] + (1.0+coef_y_b)*rhs_x[ijk] + coef_y_a*rhs_x[idx_jm]);
                }
            }
        }

        // rhs_z ---------------------------
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                for (k=1; k<nz1-1; ++k) {
                    ijk    = idx_ijk(i, j, k  , nx1, ny1);
                    idx_kp = idx_ijk(i, j, k+1, nx1, ny1);
                    idx_km = idx_ijk(i, j, k-1, nx1, ny1);
                    dzdz = DZ[k]*DZ[k];

                    coef_z_a = (dt / 2.0 / dzdz) * ( 1.0 );
                    coef_z_b = (dt / 2.0 / dzdz) * (-2.0 );
                    coef_z_c = (dt / 2.0 / dzdz) * ( 1.0 );
                    
                    rhs_z[ijk] = (coef_z_c*rhs_y[idx_kp] + (1.0+coef_z_b)*rhs_y[ijk] + coef_z_a*rhs_y[idx_km]);

                    rhs_z[ijk] += dt * ( 3.0 * Pi*Pi * cos(Pi*X[i]) * cos(Pi*Y[j]) * cos(Pi*Z[k]));
                }
            }
        }

        // Calculating A matrix ----------------------------------------------------------------

        // bdy(z) = Periodic

        // z solve
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                for (k=1; k<nz1-1; ++k) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    ki = idx_ki(k-1, i-1, nz1-2);
                    dzdz = DZ[k]*DZ[k];

                    coef_z_a = (dt / 2.0 / dzdz) * ( 1.0 );
                    coef_z_b = (dt / 2.0 / dzdz) * (-2.0 );
                    coef_z_c = (dt / 2.0 / dzdz) * ( 1.0 );

                    Azz[ki] = -coef_z_a;
                    Bzz[ki] = (1.0-coef_z_b);
                    Czz[ki] = -coef_z_c;
                    Dzz[ki] = rhs_z[ijk];
                }
            }
            tdma_cycl_many(Azz, Bzz, Czz, Dzz, nx1-2, nz1-2);
            for (i=1; i<nx1-1; ++i) {
                for (k=1; k<nz1-1; ++k) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    ki = idx_ki(k-1, i-1, nz1-2);

                    theta_z[ijk] = Dzz[ki];
                }
            }
        }

        // bdy(y)
        for (k=1; k<nz1-1; ++k) {
            for (i=1; i<nx1-1; ++i) {

                dxdx = DX[i]*DX[i];
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 );

                // j=0
                dydy = DY[0]*DY[0];
                coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) );

                idx    = idx_ik(i  , k, nx1);
                idx_ip = idx_ik(i+1, k, nx1);
                idx_im = idx_ik(i-1, k, nx1);
                theta_z[idx_ijk(i, 1, k, nx1, ny1)] += coef_y_a * (-coef_x_a*y_left_bdy[idx_im] + (1.0-coef_x_b)*y_left_bdy[idx] - coef_x_c*y_left_bdy[idx_ip]);

                // j=ny1-1
                dydy = DY[ny1-1]*DY[ny1-1];
                coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) );

                idx    = idx_ik(i  , k, nx1);
                idx_ip = idx_ik(i+1, k, nx1);
                idx_im = idx_ik(i-1, k, nx1);
                theta_z[idx_ijk(i, ny1-2, k, nx1, ny1)] += coef_y_c * (-coef_x_a*y_right_bdy[idx_im] + (1.0-coef_x_b)*y_right_bdy[idx] - coef_x_c*y_right_bdy[idx_ip]);
            }
        }

        // y solve
        for (i=1; i<nx1-1; ++i) {
            for (k=1; k<nz1-1; ++k) {
                for (j=1; j<ny1-1; ++j) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    jk = idx_jk(j-1, k-1, ny1-2);
                    dydy = DY[j]*DY[j];

                    coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * theta_y_left_index[j] + (1.0/3.0) * theta_y_right_index[j] );
                    coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -     (2.0) * theta_y_right_index[j] );
                    coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] + (5.0/3.0) * theta_y_right_index[j] );

                    Ayy[jk] = -coef_y_a;
                    Byy[jk] = (1.0-coef_y_b);
                    Cyy[jk] = -coef_y_c;
                    Dyy[jk] = theta_z[ijk];
                }
            }
            tdma_many(Ayy, Byy, Cyy, Dyy, nz1-2, ny1-2);
            for (k=1; k<nz1-1; ++k) {
                for (j=1; j<ny1-1; ++j) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    jk = idx_jk(j-1, k-1, ny1-2);

                    theta_y[ijk] = Dyy[jk];
                }
            }    
        }

        // bdy(x) = periodic

        // x solve
        for (k=1; k<nz1-1; ++k) {
            for (j=1; j<ny1-1; ++j) {
                for (i=1; i<nx1-1; ++i) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    ij = idx_ij(i-1, j-1, nx1-2);
                    dxdx = (DX[i]*DX[i]);

                    coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 );
                    coef_x_b = (dt / 2.0 / dxdx) * (-2.0 );
                    coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 );

                    Axx[ij] = -coef_x_a;
                    Bxx[ij] = (1.0-coef_x_b);
                    Cxx[ij] = -coef_x_c;
                    Dxx[ij] = theta_y[ijk];
                }
            }
            tdma_cycl_many(Axx, Bxx, Cxx, Dxx, ny1-2, nx1-2);
            for (j=1; j<ny1-1; ++j) {
                for (i=1; i<nx1-1; ++i) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    ij = idx_ij(i-1, j-1, nx1-2);

                    theta[ijk] = Dxx[ij];
                }
            }
        }
    }

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "elapsed time: " << elapsed_ms << " ms\n";

    // save results
    std::cout << "tN = " << dt * max_iter << std::endl;
    // save_3d_to_csv(theta, nx1, ny1, nz1, "results", "theta_single", 15);

    // call error
    double exact_value;
    double error = 0;

    for (int k=1; k<nz1-1; ++k) {
        for (int j=1; j<ny1-1; ++j) {
            for (int i=1; i<nx1-1; ++i) {
                exact_value = sin(Pi*X[i])*sin(Pi*Y[j])*sin(Pi*Z[k]) * exp(-3*Pi*Pi * 0.005) +
                              cos(Pi*X[i])*cos(Pi*Y[j])*cos(Pi*Z[k]);

                ijk = idx_ijk(i, j, k, nx1, ny1);
                error += pow(theta[ijk]-exact_value, 2);
            }
        }
    }
    std::cout << sqrt(error / nx1 / ny1 / nz1) << std::endl;

    return 0;
}