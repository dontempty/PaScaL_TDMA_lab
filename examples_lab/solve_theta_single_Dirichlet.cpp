// solve_theta.cpp
#include "solve_theta.hpp"
#include "../examples_lab/save.hpp"
#include "iostream"
#include <cmath>
#include <chrono> 

void solve_theta::solve_theta_plan_single(std::vector<double>& theta) 
{   
    int myrank, ierr;

    // 내가 사용할거
    int idx;
    double dxdx, dydy;
    double coef_x_a, coef_x_b, coef_x_c;
    double coef_y_a, coef_y_b, coef_y_c;
    
    // Loop and index variables
    int i, j;
    int ny1 = sub.ny_sub+1; // number of cell with ghost cell in y axis (그냥 셀 개수)
    int nx1 = sub.nx_sub+1; // number of cell with ghost cell in x axis

    // 1) 계획(plan) 객체 선언
    PaScaL_TDMA::ptdma_plan_single px_single, py_single;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    auto cx = topo.commX();
    int rankx = cx.myrank;
    PaScaL_TDMA tdma_x;
    tdma_x.PaScaL_TDMA_plan_single_create(px_single, cx.myrank, cx.nprocs, cx.comm, 0);
    std::vector<double> Axx((nx1-2)*(ny1-2)), Bxx((nx1-2)*(ny1-2)), Cxx((nx1-2)*(ny1-2)), Dxx((nx1-2)*(ny1-2));
    std::vector<double> Ax(nx1-2), Bx(nx1-2), Cx(nx1-2), Dx(nx1-2);

    auto cy = topo.commY();
    int ranky = cy.myrank;
    PaScaL_TDMA tdma_y;
    tdma_y.PaScaL_TDMA_plan_single_create(py_single, cy.myrank, cy.nprocs, cy.comm, 0);
    std::vector<double> Ayy((ny1-2)*(nx1-2)), Byy((ny1-2)*(nx1-2)), Cyy((ny1-2)*(nx1-2)), Dyy((ny1-2)*(nx1-2));
    std::vector<double> Ay(ny1-2), By(ny1-2), Cy(ny1-2), Dy(ny1-2);
    
    std::vector<double> rhs(nx1 * ny1, 0.0);
    std::vector<double> theta_half(nx1 * ny1, 0.0);

    // Calculating rhs -----------------------------------------------------------------------------------------------------

    for (j=1; j<ny1-1; ++j) {
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;
            rhs[idx] = sin(Pi*sub.x_sub[i]) * sin(Pi*sub.y_sub[j]);
        }
    }
    // save_rhs_to_csv(rhs, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 17);
    // Calculating rhs -----------------------------------------------------------------------------------------------------

    

    auto start = std::chrono::steady_clock::now();
    // Calculating A system ----------------------------------------------------------------------
    // ------------------------------------------ solve ------------------------------------------

    // x -------------------------------------------
    for (j=1; j<ny1-1; ++j) {
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;
            dxdx = (sub.dmx_sub[i]*sub.dmx_sub[i]);
            
            coef_x_a = (1.0 / dxdx) * ( 1.0 + (1.0) * sub.theta_x_right_index[i] ) * ( 1.0 - sub.theta_x_left_index[i] );
            coef_x_b = (1.0 / dxdx) * (-2.0 -     (1.0) * sub.theta_x_left_index[i] -   (1.0) * sub.theta_x_right_index[i] );
            coef_x_c = (1.0 / dxdx) * ( 1.0 + (1.0) * sub.theta_x_left_index[i] ) * ( 1.0 - sub.theta_x_right_index[i] );

            Ax[(i-1)] = coef_x_a;
            Bx[(i-1)] = coef_x_b;
            Cx[(i-1)] = coef_x_c;
            Dx[(i-1)] = rhs[idx];
        }
        tdma_x.PaScaL_TDMA_single_solve(px_single, Ax, Bx, Cx, Dx, nx1-2);
        // Return the solution to the r.h.s. line-by-line.
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;
            theta_half[idx] = Dx[i-1];
        }
    }
   
    // save_rhs_to_csv(theta_half, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 17);

    // y -------------------------------------------
    for (i=1; i<nx1-1; ++i) {
        for (j=1; j<ny1-1; ++j) {
            idx = j * nx1 + i;
            dydy = (sub.dmy_sub[j]*sub.dmy_sub[j]);
            
            coef_y_a = (1.0 / dydy) * ( 1.0 + (1.0) * sub.theta_y_right_index[j] ) * ( 1.0 - sub.theta_y_left_index[j] );
            coef_y_b = (1.0 / dydy) * (-2.0 -     (1.0) * sub.theta_y_left_index[j] -   (1.0) * sub.theta_y_right_index[j] );
            coef_y_c = (1.0 / dydy) * ( 1.0 + (1.0) * sub.theta_y_left_index[j] ) * ( 1.0 - sub.theta_y_right_index[j] );
            
            Ay[(j-1)] = coef_y_a;
            By[(j-1)] = coef_y_b;
            Cy[(j-1)] = coef_y_c;
            Dy[(j-1)] = theta_half[idx];
            
        }
        tdma_y.PaScaL_TDMA_single_solve(py_single, Ay, By, Cy, Dy, ny1-2);
        // Return the solution to the r.h.s. line-by-line.
        for (j=1; j<ny1-1; ++j) {
            idx = j * nx1 + i;
            theta[idx] = Dy[j-1];
        }  
    }
    
    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (myrank==0) {
        std::cout << "elapsed time: " << elapsed_ms << " ms\n";
    }

    tdma_x.PaScaL_TDMA_plan_single_destroy(px_single);
    tdma_y.PaScaL_TDMA_plan_single_destroy(py_single);

    // save results
    save_rhs_to_csv(theta, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 17);
}