// solve_theta.cpp
#include "solve_theta.hpp"
#include "../examples_lab/save.hpp"
#include "iostream"
#include <cmath>
#include <chrono> 

void solve_theta::solve_theta_plan_single(std::vector<double>& theta) 
{   
    int myrank, ierr;
    int idx;
    double dxdx, dydy;
    double coef_x_a, coef_x_b, coef_x_c;
    double coef_y_a, coef_y_b, coef_y_c;
    
    // Loop and index variables
    int i, j;
    int jem, jep;
    int ny1 = sub.ny_sub+1; // number of cell with ghost cell in y axis (그냥 셀 개수)
    int nx1 = sub.nx_sub+1; // number of cell with ghost cell in x axis

    // 1) 계획(plan) 객체 선언
    PaScaL_TDMA::ptdma_plan_many px_many, py_many;

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    auto cx = topo.commX();
    int rankx = cx.myrank;
    PaScaL_TDMA tdma_x;
    tdma_x.PaScaL_TDMA_plan_many_create(px_many, sub.ny_sub-1, cx.myrank, cx.nprocs, cx.comm);
    std::vector<double> Axx((nx1-2)*(ny1-2)), Bxx((nx1-2)*(ny1-2)), Cxx((nx1-2)*(ny1-2)), Dxx((nx1-2)*(ny1-2));

    auto cy = topo.commY();
    int ranky = cy.myrank;
    PaScaL_TDMA tdma_y;
    tdma_y.PaScaL_TDMA_plan_many_create(py_many, sub.nx_sub-1, cy.myrank, cy.nprocs, cy.comm);
    std::vector<double> Ayy((ny1-2)*(nx1-2)), Byy((ny1-2)*(nx1-2)), Cyy((ny1-2)*(nx1-2)), Dyy((ny1-2)*(nx1-2));
    
    std::vector<double> rhs(nx1 * ny1, 0.0);
    std::vector<double> theta_half(nx1 * ny1, 0.0);

    std::vector<double> theta_vec(nx1 * ny1);

    

    // Calculating rhs -----------------------------------------------------------------------------------------------------

    for (j=1; j<ny1-1; ++j) {
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;

            rhs[idx] = sin(Pi*sub.x_sub[i]) * sin(Pi*sub.y_sub[j]);
        }
    }

    // Calculating rhs -----------------------------------------------------------------------------------------------------

    auto start = std::chrono::steady_clock::now();

    // cal -----------------  A system 만들기 ------------------------------------------------

    // x -------------------------------------------
    for (j=1; j<ny1-1; ++j) {
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;
            dxdx = (sub.dmx_sub[i]*sub.dmx_sub[i]);
            
            coef_x_a = (1.0 / dxdx) * ( 1.0 + (1.0) * sub.theta_x_right_index[i] ) * ( 1.0 - sub.theta_x_left_index[i] );
            coef_x_b = (1.0 / dxdx) * (-2.0 -     (1.0) * sub.theta_x_left_index[i] -   (1.0) * sub.theta_x_right_index[i] );
            coef_x_c = (1.0 / dxdx) * ( 1.0 + (1.0) * sub.theta_x_left_index[i] ) * ( 1.0 - sub.theta_x_right_index[i] );

            Axx[(j-1)*(nx1-2) + (i-1)] = coef_x_a;
            Bxx[(j-1)*(nx1-2) + (i-1)] = coef_x_b;
            Cxx[(j-1)*(nx1-2) + (i-1)] = coef_x_c;
            Dxx[(j-1)*(nx1-2) + (i-1)] = rhs[idx];

            // Axx[(j-1)*(nx1-2) + (i-1)] = 1.0 * ( 1.0 - sub.theta_x_left_index[i] );
            // Bxx[(j-1)*(nx1-2) + (i-1)] = 4.0;
            // Cxx[(j-1)*(nx1-2) + (i-1)] = 1.0 * ( 1.0 - sub.theta_x_right_index[i] );
            // Dxx[(j-1)*(nx1-2) + (i-1)] = 6.0 - sub.theta_x_left_index[i] - sub.theta_x_right_index[i];
        }
    }

    tdma_x.PaScaL_TDMA_many_solve(px_many, Axx, Bxx, Cxx, Dxx, ny1-2, nx1-2);
    for (j=1; j<ny1-1; ++j) {
        for (i=1; i<nx1-1; ++i) {
            idx = j * nx1 + i;
            theta_half[idx] = Dxx[(j-1)*(nx1-2) + (i-1)]; 
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

            Ayy[(i-1)*(ny1-2) + (j-1)] = coef_y_a;
            Byy[(i-1)*(ny1-2) + (j-1)] = coef_y_b;
            Cyy[(i-1)*(ny1-2) + (j-1)] = coef_y_c;
            Dyy[(i-1)*(ny1-2) + (j-1)] = theta_half[idx];
        }          
    }
    tdma_y.PaScaL_TDMA_many_solve(py_many, Ayy, Byy, Cyy, Dyy, nx1-2, ny1-2);
    for (i=1; i<nx1-1; ++i) {
        for (j=1; j<ny1-1; ++j) {
            idx = j * nx1 + i;
            theta[idx] = Dyy[(i-1)*(ny1-2) + (j-1)];
        }
    }

    // cal -----------------  A system 만들기 ------------------------------------------------

    auto end = std::chrono::steady_clock::now();
    auto elapsed_ms = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count();
    if (myrank==0) {
        std::cout << "elapsed time: " << elapsed_ms << " ms\n";
    }

    tdma_x.PaScaL_TDMA_plan_many_destroy(px_many, px_many.nprocs);
    tdma_y.PaScaL_TDMA_plan_many_destroy(py_many, py_many.nprocs);

    // save results
    save_rhs_to_csv(theta, nx1, ny1, "results", "rhs_" + std::to_string(cy.myrank) + std::to_string(cx.myrank) +".csv", 17);
}