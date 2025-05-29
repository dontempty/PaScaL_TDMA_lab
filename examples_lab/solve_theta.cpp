// solve_theta.cpp
#include "solve_theta.hpp"
#include "iostream"

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
    
    // 일단 포아송 방정식 푼다. 
    // u = sin(pi*x) * sin(pi*y)
    // f = -2*pi * sin(pi*x) * sin(pi*y)
    
    // Calculating r.h.s -----------------------------------------------------------------------------------------------------
    // std::vector<double> rhs((sub.ny_sub+1) * (sub.nx_sub));

    // y axis ----------------------------------------------------------------------------------------------------------------

    // PaScaL_TDMA tdma_y;
    // auto cy = topo.commY();
    // tdma_x.PaScaL_TDMA_plan_single_create(px_single, cy.myrank, cy.nprocs, cy.comm, 0);
    
    // x axis ----------------------------------------------------------------------------------------------------------------
    // PaScaL_TDMA tdma_x;
    // auto cx = topo.commX();
    
    // tdma_x.PaScaL_TDMA_plan_single_create(px_single, cx.myrank, cx.nprocs, cx.comm, 0);


    

}