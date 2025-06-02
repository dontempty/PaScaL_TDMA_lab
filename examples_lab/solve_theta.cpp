// solve_theta.cpp
#include "solve_theta.hpp"
#include "../examples_lab/save.hpp"
#include "iostream"
#include <cmath>

void solve_theta::solve_theta_plan_single(double* theta) 
{   
    int myrank, ierr;
    int time_step;      // Current time step
    double t_curr;      // Current simulation time

    // 내가 사용할거
    int idx;
    int idx_ip, idx_im;
    int idx_jp, idx_jm;
    double dxdx, dydy;
    double coef_x_a, coef_x_b, coef_x_c;
    double coef_y_a, coef_y_b, coef_y_c;
    
    // Loop and index variables
    int i, j;
    int ip, jp;
    int im, jm;
    int jem, jep;
    int ny1 = sub.ny_sub+1; // number of cell with ghost cell in y axis (그냥 셀 개수)
    int nx1 = sub.nx_sub+1; // number of cell with ghost cell in x axis

    // Temporary variables for coefficient computations
    double dedx1, dedx2, dedy3, dedy4, dedz5, dedz6;     // Derivative terms
    double viscous_e1, viscous_e2, viscous_e3, viscous;  // Viscous terms
    double ebc_down, ebc_up, ebc;                        // Boundary terms
    double eAPI, eAMI, eACI;                             // Diffusion treatment terms in x-direction
    double eAPJ, eAMJ, eACJ;                             // Diffusion treatment terms in y-direction
    double eAPK, eAMK, eACK;                             // Diffusion treatment terms in z-direction
    double eRHS;                                         // From eAPI to eACK

    // 1) 계획(plan) 객체 선언
    PaScaL_TDMA::ptdma_plan_single px_single, py_single;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    // t_curr = tStart 나중에 heat eq 풀면 그 때 사용하자
    // dt = dtstart

    if (myrank==0) {
        std::cout << "Start to solve" << std::endl;
    }

    // 4. [ ][ ] | [ ][ ] | [ ][ ]
    // 3. [ ][ ] | [ ][ ] | [ ][ ]
    //    ------------------------
    // 2. [ ][ ] | [ ][ ] | [ ][ ]
    // 1. [ ][ ] | [ ][ ] | [ ][ ]    일단 도메인을 이렇게 분해했다고 가정

    // 우선 1. 번 기준으로 rhs, A 만들고 single로 풀기
    // 그런 다음 위 과정을 2번에서 반복


    // 일단 포아송 방정식 푼다. 
    // BCFD - HW7 을 푼다.

    auto cx = topo.commX();
    int rankx = cx.myrank;
    PaScaL_TDMA tdma_x;
    tdma_x.PaScaL_TDMA_plan_single_create(px_single, cx.myrank, cx.nprocs, cx.comm, 0);
    std::vector<double> Ax(nx1-2), Bx(nx1-2), Cx(nx1-2), Dx(nx1-2);

    auto cy = topo.commY();
    int ranky = cy.myrank;
    PaScaL_TDMA tdma_y;
    tdma_y.PaScaL_TDMA_plan_single_create(py_single, cy.myrank, cy.nprocs, cy.comm, 0);
    std::vector<double> Ay(ny1-2), By(ny1-2), Cy(ny1-2), Dy(ny1-2);
    
    std::vector<double> rhs_x(nx1 * ny1, 0.0);
    std::vector<double> rhs_y(nx1 * ny1, 0.0);
    std::vector<double> theta_old(nx1 * ny1, 0.0);
    
    double dt = 0.01;
    int max_iter = 1000;
    int check_num = 25;
    double tol = 1e-12;
    double error = 0;
    double global_error = 0.0;

    for (int t_step=0; t_step<max_iter; ++t_step) {

        // ---------------------- 차분 공식 -----------------------------
        // bdy 에서 거리 정보가 달라진다.
        // ( u(i+1) - 3u(i) + 2u(i-1) ) / (3/4)dxdx  left_bdy
        // ( u(i+1) - 2u(i) + u(i-1) ) / (3/4)dxdx   interior
        // ( 2u(i+1) - 3u(i) + u(i-1) ) / (3/4)dxdx  right_bdy
        // -----------------------------

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

        // (1+Ax/2)(1+Ay/2)u_n
        // y ---------
        for (j=1; j<ny1-1; ++j) {
            for (i=0; i<nx1; ++i) {
                idx = j * nx1 + i;
                idx_jm = (j-1) * nx1 + i;
                idx_jp = (j+1) * nx1 + i;
                dydy = (sub.dmy_sub[j]*sub.dmy_sub[j]);

                coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (5.0/3.0) * sub.theta_y_left_index[j] + (1.0/3.0) * sub.theta_y_right_index[j] );
                coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * sub.theta_y_left_index[j] -     (2.0) * sub.theta_y_right_index[j] );
                coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * sub.theta_y_left_index[j] + (5.0/3.0) * sub.theta_y_right_index[j] );

                rhs_y[idx] += (coef_y_c*theta[idx_jp] + (1+coef_y_b)*theta[idx] + coef_y_a*theta[idx_jm]);
            }
        }
        
        // x ---------
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                idx_im = j * nx1 + (i-1);
                idx_ip = j * nx1 + (i+1);
                dxdx = (sub.dmx_sub[i]*sub.dmx_sub[i]);
                
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (5.0/3.0) * sub.theta_x_left_index[i] + (1.0/3.0) * sub.theta_x_right_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * sub.theta_x_left_index[i] -     (2.0) * sub.theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * sub.theta_x_left_index[i] + (5.0/3.0) * sub.theta_x_right_index[i] );

                rhs_x[idx] += (coef_x_c*rhs_y[idx_ip] + (1+coef_x_b)*rhs_y[idx] + coef_x_a*rhs_y[idx_im]);
                
                // source func (S = 2(2-x^2-y^2))
                rhs_x[idx] += (dt) * 2.0 * (2.0 - sub.x_sub[i]*sub.x_sub[i] - sub.y_sub[j]*sub.y_sub[j]);
                // rhs_x[idx] += dt * sin(Pi * sub.x_sub[i]) * sin(Pi * sub.y_sub[j]);
            }
        }

        // cal -----------------  A system 만들기 ------------------------------------------------

        // x -------------------------------------------
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                dxdx = (sub.dmx_sub[i]*sub.dmx_sub[i]);
                
                coef_x_a = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * sub.theta_x_right_index[i] ) * ( 1.0 - sub.theta_x_left_index[i] );
                coef_x_b = (dt / 2.0 / dxdx) * (-2.0 -     (2.0) * sub.theta_x_left_index[i] -   (2.0) * sub.theta_x_right_index[i] );
                coef_x_c = (dt / 2.0 / dxdx) * ( 1.0 + (1.0/3.0) * sub.theta_x_left_index[i] ) * ( 1.0 - sub.theta_x_right_index[i] );

                Ax[i-1] = -coef_x_a;
                Bx[i-1] = 1-coef_x_b;
                Cx[i-1] = -coef_x_c;
                Dx[i-1] = rhs_x[idx];
            }
            tdma_x.PaScaL_TDMA_single_solve(px_single, Ax, Bx, Cx, Dx, nx1-2);
            // Return the solution to the r.h.s. line-by-line.
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                rhs_x[idx] = Dx[i-1];
            }
        }


        // y -------------------------------------------
        for (i=1; i<nx1-1; ++i) {
            for (j=1; j<ny1-1; ++j) {
                idx = j * nx1 + i;
                dydy = (sub.dmy_sub[j]*sub.dmy_sub[j]);
                
                coef_y_a = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * sub.theta_y_right_index[j] ) * ( 1.0 - sub.theta_y_left_index[j] );
                coef_y_b = (dt / 2.0 / dydy) * (-2.0 -     (2.0) * sub.theta_y_left_index[j] -   (2.0) * sub.theta_y_right_index[j] );
                coef_y_c = (dt / 2.0 / dydy) * ( 1.0 + (1.0/3.0) * sub.theta_y_left_index[j] ) * ( 1.0 - sub.theta_y_right_index[j] );

                Ay[j-1] = -coef_y_a;
                By[j-1] = 1-coef_y_b;
                Cy[j-1] = -coef_y_c;
                Dy[j-1] = rhs_x[idx];
            }
            tdma_y.PaScaL_TDMA_single_solve(py_single, Ay, By, Cy, Dy, ny1-2);
            // Return the solution to the theta. line-by-line.
            for (j=1; j<ny1-1; ++j) {
                idx = j * nx1 + i;
                theta[idx] = Dy[j-1];
            }
        }

        // Update ghostcells from the solutions.
        sub.ghostcellUpdate(theta, cx, cy, params);

        // break point ----------------------
        error = 0;
        for (j=1; j<ny1-1; ++j) {
            for (i=1; i<nx1-1; ++i) {
                idx = j * nx1 + i;
                error += (theta_old[idx]-theta[idx])*(theta_old[idx]-theta[idx]);
            }
        }
        error  = sqrt ( error / ( (ny1-2)*(nx1-2) ) );

        MPI_Allreduce(&error, &global_error, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (myrank==0 && t_step%check_num == 0) {
            std::cout << "Step = " << t_step << "| global_error = " << global_error << std::endl;
        }

        if (global_error < tol) {
            std::cout << "Converge in " << t_step << "step" << std::endl;
            break;
        }
   
    }   // Time step end------------------------
    tdma_x.PaScaL_TDMA_plan_single_destroy(px_single);
    tdma_y.PaScaL_TDMA_plan_single_destroy(py_single);

    // save results
    std::vector<double> theta_vec(nx1 * ny1);
    for (int j = 0; j < ny1; ++j) {
        for (int i = 0; i < nx1; ++i) {
            theta_vec[j*nx1 + i] = theta[j*nx1 + i];
        }
    }
    save_rhs_to_csv(theta_vec, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 13);

    // ---------------- debug ------------------------------------------------
    // std::vector<double> rhs(nx1 * ny1, 0.0);
    // for (j=1; j<ny1-1; ++j) {
    //     for (i=1; i<nx1-1; ++i) {
    //         idx = j * nx1 + i;
    //         if (cx.myrank==0 && i==1) {
    //             rhs[idx] = 5;
    //         } else if (cx.myrank==1 && i==nx1-2) {
    //             rhs[idx] = 5;
    //         } else {

    //             rhs[idx] = 6;
    //         }
    //     }
    // }
    
    // for (j=1; j<ny1-1; ++j) {
    //     for (i=1; i<nx1-1; ++i) {
    //         idx = j * nx1 + i;

    //         Ax[i-1] = 1;
    //         Bx[i-1] = 4;
    //         Cx[i-1] = 1;
    //         Dx[i-1] = rhs[idx];
    //     }
    //     tdma_x.PaScaL_TDMA_single_solve(px_single, Ax, Bx, Cx, Dx, nx1-2);
    //     for (i=1; i<nx1-1; ++i) {
    //         idx = j * nx1 + i;
    //         theta[idx] = Dx[i-1];
    //     }
    // }

    // std::vector<double> theta_vec(nx1 * ny1);
    // for (int j = 0; j < ny1; ++j) {
    //     for (int i = 0; i < nx1; ++i) {
    //         theta_vec[j*nx1 + i] = theta[j*nx1 + i];
    //     }
    // }
    // save_rhs_to_csv(theta_vec, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 13);

    // 방금 이걸로 디버깅 했는데 존나 잘 푼다.
    // 뭐가 문제인거임 ???????????????????????????????????
}