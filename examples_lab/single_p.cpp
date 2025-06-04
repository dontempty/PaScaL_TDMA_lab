#include "save.hpp"
#include "iostream"
#include <vector>
#include <cmath>
#include <chrono> 

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

int main() {
    int i, j;

    int Nx = 1024; 
    int Ny = 1024;
    int nx1 = Nx+2;
    int ny1 = Ny+2;

    double x0 = -1;
    double xN = 1;
    double y0 = -1;
    double yN = 1;

    double dy = (yN-y0)/Ny;
    double dx = (xN-x0)/Nx;

    std::vector<double> X(nx1);
    std::vector<double> Y(ny1);
    for (i=0; i<nx1; ++i) {
        if (i==0) {
            X[i] = x0;
            Y[i] = y0;
        }
        else if (i==nx1-1) {
            X[i] = xN;
            Y[i] = yN;
        }
        else {
            X[i] = x0 + dx/2 + (i-1)*dx;
            Y[i] = y0 + dy/2 + (i-1)*dy;
        }
    }

    std::vector<double> theta_x_left_index(nx1);
    std::vector<double> theta_x_right_index(nx1);
    std::vector<double> theta_y_left_index(ny1);
    std::vector<double> theta_y_right_index(ny1);
    theta_x_left_index[1] = 1;
    theta_x_right_index[nx1-2] = 1;
    theta_y_left_index[1] = 1;
    theta_y_right_index[nx1-2] = 1;

    std::vector<double> Ax(nx1-2), Bx(nx1-2), Cx(nx1-2), Dx(nx1-2);
    std::vector<double> Ay(ny1-2), By(ny1-2), Cy(ny1-2), Dy(ny1-2);

    std::vector<double> rhs_x(nx1 * ny1, 0.0);
    std::vector<double> rhs_y(nx1 * ny1, 0.0);
    std::vector<double> theta(nx1 * ny1, 0.0);
    std::vector<double> theta_half(nx1 * ny1, 0.0);
    std::vector<double> theta_old(nx1 * ny1, 0.0);

    std::vector<double> theta_vec(nx1 * ny1);
    

    int idx;
    int idx_ip, idx_im;
    int idx_jp, idx_jm;
    double dxdx, dydy;
    double coef_x_a, coef_x_b, coef_x_c;
    double coef_y_a, coef_y_b, coef_y_c;

    double dt = 0.01;
    int max_iter = 1000;
    int check_num = 50;
    double tol = 1e-12;
    double error = 0;
    double global_error = 0.0;

    auto start = std::chrono::steady_clock::now();
    for (int t_step=0; t_step<max_iter; ++t_step) {

        // copy theta to theta_old
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                theta_old[idx] = theta[idx];
            }
        }

        // Calculating r.h.s -----------------------------------------------------------------------------------------------------

        // rhs init
        std::fill(rhs_x.begin(), rhs_x.end(), 0.0);
        std::fill(rhs_y.begin(), rhs_y.end(), 0.0);

        // y ---------
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                idx_jm = (j-1) * nx1 + i;
                idx_jp = (j+1) * nx1 + i;
                dydy = dy*dy;

                coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * theta_y_left_index[j] + (1.0/3.0) * theta_y_right_index[j] );
                coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -     (2.0) * theta_y_right_index[j] );
                coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] + (5.0/3.0) * theta_y_right_index[j] );

                rhs_y[idx] += (coef_y_c*theta[idx_jp] + (1+coef_y_b)*theta[idx] + coef_y_a*theta[idx_jm]);
            }
        }

        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                idx_im = j * nx1 + (i-1);
                idx_ip = j * nx1 + (i+1);
                dxdx = dx*dx;
                
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * theta_x_left_index[i] + (1.0/3.0) * theta_x_right_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -     (2.0) * theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] + (5.0/3.0) * theta_x_right_index[i] );

                rhs_x[idx] += (coef_x_c*rhs_y[idx_ip] + (1+coef_x_b)*rhs_y[idx] + coef_x_a*rhs_y[idx_im]);
                
                // source func (S = 2(2-x^2-y^2))
                rhs_x[idx] += (dt) * 2.0 * (2.0 - X[i]*X[i] - Y[j]*Y[j]);
            }
        }

        //  solve ----------
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                dxdx = dx*dx;
                
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_right_index[i] ) * ( 1.0 - theta_x_left_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * theta_x_left_index[i] -   (2.0) * theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * theta_x_left_index[i] ) * ( 1.0 - theta_x_right_index[i] );

                Ax[i-1] = -coef_x_a;
                Bx[i-1] = 1-coef_x_b;
                Cx[i-1] = -coef_x_c;
                Dx[i-1] = rhs_x[idx];
            }
            tdma_single(Ax, Bx, Cx, Dx, nx1-2);
            // Return the solution to the r.h.s. line-by-line.
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                theta_half[idx] = Dx[i-1];
            }
        }

        for (i=1; i<nx1-1; ++i) {
            for (j=1; j<ny1-1; ++j) {
                idx = j * nx1 + i;
                dydy = dy*dy;
                
                coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_right_index[j] ) * ( 1.0 - theta_y_left_index[j] );
                coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * theta_y_left_index[j] -   (2.0) * theta_y_right_index[j] );
                coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * theta_y_left_index[j] ) * ( 1.0 - theta_y_right_index[j] );

                Ay[j-1] = -coef_y_a;
                By[j-1] = 1-coef_y_b;
                Cy[j-1] = -coef_y_c;
                Dy[j-1] = theta_half[idx];
            }
            tdma_single(Ay, By, Cy, Dy, ny1-2);
            // Return the solution to the theta. line-by-line.
            for (j=1; j<ny1-1; ++j) {
                idx = j * nx1 + i;
                theta[idx] = Dy[j-1];
            }            
        }

        // break point ----------------------
        error = 0;
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                error += (theta_old[idx]-theta[idx])*(theta_old[idx]-theta[idx]);
            }
        }
        global_error  = sqrt ( error / ( (ny1-2)*(nx1-2) ) );

        if (t_step%check_num == 0) {
            std::cout << "Step = " << t_step << "| global_error = " << global_error << std::endl;
        }

        if (global_error < tol) {
            std::cout << "Converge in " << t_step << "step" << std::endl;
            break;
        }
    }
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    std::cout << "elapsed time: " << elapsed_ms << " ms\n";

    for (int j = 0; j < ny1; ++j) {
        for (int i = 0; i < nx1; ++i) {
            theta_vec[j*nx1 + i] = theta[j*nx1 + i];
        }
    }
    save_rhs_to_csv(theta_vec, nx1, ny1, "results", "rhs_single.csv", 13);

    return 0;
}