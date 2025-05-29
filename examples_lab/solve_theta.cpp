// solve_theta.cpp
#include "solve_theta.hpp"
#include "iostream"
#include <cmath>

void solve_theta::solve_theta_plan_single(double* theta) 
{   
    int myrank, ierr;
    int time_step;      // Current time step
    double t_curr;      // Current simulation time
    
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
    // u = sin(pi*x) * sin(pi*y)
    // f = -2*pi * sin(pi*x) * sin(pi*y)

    int idx;

    // Calculating r.h.s -----------------------------------------------------------------------------------------------------
    std::vector<double> rhs(nx1 * ny1);
    for (i=0; i<nx1; ++i) {
        for (j=0; j<ny1; ++j) {
            idx = i * ny1 + j;
            rhs[idx] = -2*Pi*sin(Pi*sub.x_sub[i])*sin(Pi*sub.y_sub[j]);     // source func
            rhs[idx] += params.theta_x_L_D * sub.theta_x_left_index[i];     // Drichlet bdy  x_left
            rhs[idx] += params.theta_x_R_D * sub.theta_x_right_index[i];    // Drichlet bdy  x_right
            rhs[idx] += params.theta_y_L_D * sub.theta_y_left_index[j];     // Drichlet bdy  y_left
            rhs[idx] += params.theta_y_R_D * sub.theta_y_right_index[j];    // Drichlet bdy  y_right
        }
    }
    // std::cout << "myrank: " << myrank << "| n_sub(xy)" << nx1 << ny1 << "|";
    // for (int i=0; i<nx1*ny1; ++i) {
    //     std::cout << rhs[i] << " ";
    // }

    // ---------------------- 차분 공식 -----------------------------
    // ( u(i+1) - 3u(i) + 2u(i-1) ) / dxdx  left_bdy
    // ( u(i+1) - 2u(i) + u(i-1) ) / dxdx   interior
    // ( 2u(i+1) - 3u(i) + u(i-1) ) / dxdx  right_bdy
    // -----------------------------


    // x axis ----------------------------------------------------------------------------------------------------------------
    auto cx = topo.commX();
    PaScaL_TDMA tdma_x;
    tdma_x.PaScaL_TDMA_plan_single_create(px_single, cx.myrank, cx.nprocs, cx.comm, 0);

    std::vector<double> Ax(nx1-2), Bx(nx1-2), Cx(nx1-2), Dx(nx1-2);
    double dxdx;
    // for (j=1; j<ny1-1; ++j)
    for (j=1; j<2; ++j) {   // 여기는 y 층으로 올라가는 loop
        for (i=1; i<nx1-1; ++i) {   // 여기는 system 만드는 곳
            idx = i * ny1 + j;

            dxdx = (sub.dmx_sub[i]*sub.dmx_sub[i]);     // 지금은 그냥 (l/(nx-1) 로 일정하다)
            Ax[i-1] = 1 / dxdx;
            Ax[i-1] += (1 / dxdx) * sub.theta_x_left_index[i];   // A는 왼쪽인 경우만 바뀐다.

            Bx[i-1] = -2 / dxdx;
            Bx[i-1] += -(1 / dxdx) * sub.theta_x_left_index[i];
            Bx[i-1] += -(1 / dxdx) * sub.theta_x_right_index[i];

            Cx[i-1] = 1 / dxdx;
            Cx[i-1] += (1 / dxdx) * sub.theta_x_right_index[i];

            Dx[i-1] = rhs[idx];
        }
        std::cout << "myrank: " << myrank << "| n(xy)1 | " << nx1 << "," << ny1 << "|";
        for (int i=0; i<nx1-2; ++i) {
            std::cout << Cx[0] << " ";
        }

        // solve !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        // tdma_x.PaScaL_TDMA_single_solve(px_single, Ax, Bx, Cx, Dx);

    }
    

    // y axis ----------------------------------------------------------------------------------------------------------------

    // PaScaL_TDMA tdma_y;
    // auto cy = topo.commY();
    // tdma_y.PaScaL_TDMA_plan_single_create(px_single, cy.myrank, cy.nprocs, cy.comm, 0);

    // for (i=1; i)
    
    


    

}