#include "save.hpp"
#include "iostream"
#include <vector>
#include <cmath>
#include <chrono> 

#include "index.hpp"

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
        rr = 1.0 / (b[i] - a[i]*c[i-1]);
        d[i] = rr*(d[i] - a[i]*d[i-1]);
        e[i] = rr*(e[i] - a[i]*e[i-1]);
        c[i] = rr*c[i];
    }

    for (i=n1-2; i>=1; --i) {
        d[i] = d[i] - c[i]*d[i+1];
        e[i] = e[i] - c[i]*e[i+1];
    }

    d[0] = (d[0] - a[0]*d[n1-1] - c[0]*d[1])/(b[0] + a[0]*e[n1-1] + c[0]*e[1]);

    for (i=1; i<=n1-1; ++i) {
        d[i] = d[i] + d[0]*e[i];
    }
}

// compile = g++ -O2 -std=c++17 -static-libgcc -static-libstdc++ -o single_fixed_D_3D.exe single_fixed_D_3D.cpp
// Run = ./single_fixed_D_3D.exe 
int main() {

    int i, j, k;

    // paramater
    int Nx = 64; 
    int Ny = 64;
    int Nz = 64;
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

    std::vector<double> Ax(nx1-2), Bx(nx1-2), Cx(nx1-2), Dx(nx1-2);
    std::vector<double> Ay(ny1-2), By(ny1-2), Cy(ny1-2), Dy(ny1-2);
    std::vector<double> Az(nz1-2), Bz(nz1-2), Cz(nz1-2), Dz(nz1-2);

    std::vector<double> rhs_x(nx1 * ny1 * nz1, 0.0);
    std::vector<double> rhs_y(nx1 * ny1 * nz1, 0.0);
    std::vector<double> rhs_z(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta_z(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta_y(nx1 * ny1 * nz1, 0.0);
    std::vector<double> theta(nx1 * ny1 * nz1, 0.0);
    
    int idx;
    int ij, ik, jk;
    int ijk;
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

    // bdy variable
    std::vector<double> x_left_bdy(ny1 * nz1, 0.0);
    std::vector<double> x_right_bdy(ny1 * nz1, 0.0);
    for (k=0; k<nz1; ++k) {
        for (j=0; j<ny1; ++j) {
            jk = idx_jk(j, k, ny1);

            ijk = idx_ijk(0, j, k, nx1, ny1);
            x_left_bdy[jk] = theta[ijk];

            ijk = idx_ijk(nx1-1, j, k, nx1, ny1);
            x_right_bdy[jk] = theta[ijk];
        }
    }

    std::vector<double> y_left_bdy(nx1 * nz1, 0.0);
    std::vector<double> y_right_bdy(nx1 * nz1, 0.0);
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

            ijk = idx_ijk(i, j, 0, nx1, ny1);
            z_left_bdy[ij] = theta[ijk];

            ijk = idx_ijk(i, j, nz1-1, nx1, ny1);
            z_right_bdy[ij] = theta[ijk];
        }
    }

    int max_iter = 100;
    double dt = 0.001;
    int time;
    auto start = std::chrono::steady_clock::now();
    for (time=0; time<max_iter; ++time) {

        // Calculating r.h.s -------------------------------------------------------------------

        // rhs_x ---------------------------
        for (k=0; k<nz1; ++k) {
            for (j=0; j<ny1; ++j) {
                for (i=1; i<nx1-1; ++i) {
                    ijk    = idx_ijk(i  , j, k, nx1, ny1);
                    idx_ip = idx_ijk(i+1, j, k, nx1, ny1);
                    idx_im = idx_ijk(i-1, j, k, nx1, ny1);
                    dxdx = DX[i]*DX[i];

                    coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * theta_x_left_index[i] + (1.0/3.0) * theta_x_right_index[i] );
                    coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -     (2.0) * theta_x_right_index[i] );
                    coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] + (5.0/3.0) * theta_x_right_index[i] );
                    
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

                    coef_z_a = (dt / 2.0 / dzdz) * ( 1.0 + (5.0/3.0) * theta_z_left_index[k] + (1.0/3.0) * theta_z_right_index[k] );
                    coef_z_b = (dt / 2.0 / dzdz) * (-2.0 -     (2.0) * theta_z_left_index[k] -     (2.0) * theta_z_right_index[k] );
                    coef_z_c = (dt / 2.0 / dzdz) * ( 1.0 + (1.0/3.0) * theta_z_left_index[k] + (5.0/3.0) * theta_z_right_index[k] );
                    
                    rhs_z[ijk] = (coef_z_c*rhs_y[idx_kp] + (1.0+coef_z_b)*rhs_y[ijk] + coef_z_a*rhs_y[idx_km]);

                    rhs_z[ijk] += dt * ( 3.0 * Pi*Pi * cos(Pi*X[i]) * cos(Pi*Y[j]) * cos(Pi*Z[k]));
                }
            }
        }

        // Calculating A matrix ----------------------------------------------------------------

        // bdy(z)
        for (j=1; j<ny1-1; ++j) {

            dydy = DY[j]*DY[j];
            coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * theta_y_left_index[j] + (1.0/3.0) * theta_y_right_index[j] );
            coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -     (2.0) * theta_y_right_index[j] );
            coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] + (5.0/3.0) * theta_y_right_index[j] ); 

            for (i=1; i<nx1-1; ++i) {

                dxdx = DX[i]*DX[i];
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * theta_x_left_index[i] + (1.0/3.0) * theta_x_right_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -     (2.0) * theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] + (5.0/3.0) * theta_x_right_index[i] );

                // k=0
                dzdz = DZ[0]*DZ[0];
                coef_z_a = (dt / 2.0 / dzdz) * ( 1.0 + (5.0/3.0) );

                idx    = idx_ij(i  , j-1, nx1);
                idx_ip = idx_ij(i+1, j-1, nx1);
                idx_im = idx_ij(i-1, j-1, nx1);
                rhs_z[idx_ijk(i, j, 1, nx1, ny1)] += (coef_z_a) * (-coef_y_a)    * (-coef_x_a*z_left_bdy[idx_im] + (1.0-coef_x_b)*z_left_bdy[idx] - coef_x_c*z_left_bdy[idx_ip]);

                idx    = idx_ij(i  , j, nx1);
                idx_ip = idx_ij(i+1, j, nx1);
                idx_im = idx_ij(i-1, j, nx1);
                rhs_z[idx_ijk(i, j, 1, nx1, ny1)] += (coef_z_a) * (1.0-coef_y_b) * (-coef_x_a*z_left_bdy[idx_im] + (1.0-coef_x_b)*z_left_bdy[idx] - coef_x_c*z_left_bdy[idx_ip]);

                idx    = idx_ij(i  , j+1, nx1);
                idx_ip = idx_ij(i+1, j+1, nx1);
                idx_im = idx_ij(i-1, j+1, nx1);
                rhs_z[idx_ijk(i, j, 1, nx1, ny1)] += (coef_z_a) * (-coef_y_c)    * (-coef_x_a*z_left_bdy[idx_im] + (1.0-coef_x_b)*z_left_bdy[idx] - coef_x_c*z_left_bdy[idx_ip]);

                // k=nz1-1
                dzdz = DZ[nz1-1]*DZ[nz1-1];
                coef_z_c = (dt / 2.0 / dzdz) * ( 1.0 + (5.0/3.0) );

                idx    = idx_ij(i  , j-1, nx1);
                idx_ip = idx_ij(i+1, j-1, nx1);
                idx_im = idx_ij(i-1, j-1, nx1);
                rhs_z[idx_ijk(i, j, nz1-2, nx1, ny1)] += (coef_z_c) * (-coef_y_a)    * (-coef_x_a*z_right_bdy[idx_im] + (1.0-coef_x_b)*z_right_bdy[idx] - coef_x_c*z_right_bdy[idx_ip]);

                idx    = idx_ij(i  , j, nx1);
                idx_ip = idx_ij(i+1, j, nx1);
                idx_im = idx_ij(i-1, j, nx1);
                rhs_z[idx_ijk(i, j, nz1-2, nx1, ny1)] += (coef_z_c) * (1.0-coef_y_b) * (-coef_x_a*z_right_bdy[idx_im] + (1.0-coef_x_b)*z_right_bdy[idx] - coef_x_c*z_right_bdy[idx_ip]);

                idx    = idx_ij(i  , j+1, nx1);
                idx_ip = idx_ij(i+1, j+1, nx1);
                idx_im = idx_ij(i-1, j+1, nx1);
                rhs_z[idx_ijk(i, j, nz1-2, nx1, ny1)] += (coef_z_c) * (-coef_y_c)    * (-coef_x_a*z_right_bdy[idx_im] + (1.0-coef_x_b)*z_right_bdy[idx] - coef_x_c*z_right_bdy[idx_ip]);
            }
        }

        // z solve
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                for (k=1; k<nz1-1; ++k) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    dzdz = DZ[k]*DZ[k];

                    coef_z_a = (dt / 2.0 / dzdz) * ( 1.0 + (5.0/3.0) * theta_z_left_index[k] + (1.0/3.0) * theta_z_right_index[k] );
                    coef_z_b = (dt / 2.0 / dzdz) * (-2.0 -     (2.0) * theta_z_left_index[k] -     (2.0) * theta_z_right_index[k] );
                    coef_z_c = (dt / 2.0 / dzdz) * ( 1.0 + (1.0/3.0) * theta_z_left_index[k] + (5.0/3.0) * theta_z_right_index[k] );

                    Az[k-1] = -coef_z_a;
                    Bz[k-1] = (1.0-coef_z_b);
                    Cz[k-1] = -coef_z_c;
                    Dz[k-1] = rhs_z[ijk];
                }
                tdma_single(Az, Bz, Cz, Dz, nz1-2);
                for (k=1; k<nz1-1; ++k) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    theta_z[ijk] = Dz[k-1];
                }
            }
        }

        // bdy(y)
        for (k=1; k<nz1-1; ++k) {
            for (i=1; i<nx1-1; ++i) {

                dxdx = DX[i]*DX[i];
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * theta_x_left_index[i] + (1.0/3.0) * theta_x_right_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -     (2.0) * theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] + (5.0/3.0) * theta_x_right_index[i] );

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
                    dydy = DY[j]*DY[j];

                    coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * theta_y_left_index[j] + (1.0/3.0) * theta_y_right_index[j] );
                    coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -     (2.0) * theta_y_right_index[j] );
                    coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] + (5.0/3.0) * theta_y_right_index[j] );

                    Ay[j-1] = -coef_y_a;
                    By[j-1] = (1.0-coef_y_b);
                    Cy[j-1] = -coef_y_c;
                    Dy[j-1] = theta_z[ijk];
                }
                tdma_single(Ay, By, Cy, Dy, ny1-2);
                for (j=1; j<ny1-1; ++j) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    theta_y[ijk] = Dy[j-1];
                }
            }
        }

        // bdy(x)
        for (k=1; k<nz1-1; ++k) {
            for (j=1; j<ny1-1; ++j) {
                
                // i=0
                dxdx = DX[0]*DX[0];
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) );

                idx = idx_jk(j, k, ny1);
                theta_y[idx_ijk(1, j, k, nx1, ny1)] += coef_x_a * x_left_bdy[idx];

                // i=nx1-1
                dxdx = DX[nx1-1]*DX[nx1-1];
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) );

                idx = idx_jk(j, k, ny1);
                theta_y[idx_ijk(nx1-2, j, k, nx1, ny1)] += coef_x_c * x_right_bdy[idx];

            }
        }

        // x solve
        for (k=1; k<nz1-1; ++k) {
            for (j=1; j<ny1-1; ++j) {
                for (i=1; i<nx1-1; ++i) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    dxdx = (DX[i]*DX[i]);

                    coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * theta_x_left_index[i] + (1.0/3.0) * theta_x_right_index[i] );
                    coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -     (2.0) * theta_x_right_index[i] );
                    coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] + (5.0/3.0) * theta_x_right_index[i] );

                    Ax[i-1] = -coef_x_a;
                    Bx[i-1] = (1.0-coef_x_b);
                    Cx[i-1] = -coef_x_c;
                    Dx[i-1] = theta_y[ijk];
                }
                tdma_single(Ax, Bx, Cx, Dx, nx1-2);
                for (i=1; i<nx1-1; ++i) {
                    ijk = idx_ijk(i, j, k, nx1, ny1);
                    theta[ijk] = Dx[i-1];
                }
            }
        }
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "elapsed time: " << elapsed_ms << " ms\n";

    // save results
    std::cout << "tN = " << dt * max_iter << std::endl;
    save_3d_to_csv(theta, nx1, ny1, nz1, "results", "theta_single", 15);

    return 0;
}