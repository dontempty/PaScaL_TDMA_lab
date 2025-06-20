#include <mpi.h>
#include <vector>
#include "iostream"
#include "pascal_tdma.hpp"
#include "tdmas.hpp"
#include "../examples_lab/mpi_subdomain.hpp"
#include "../examples_lab/mpi_topology.hpp"

// Create a plan for a single tridiagonal system of equations.
void PaScaL_TDMA::PaScaL_TDMA_plan_single_create(ptdma_plan_single& plan, int myrank, int nprocs, MPI_Comm mpi_world, int gather_rank) {

    int nr_rd;      // Number of rows of a reduced tridiagonal system per process, 2
    int nr_rt;      // Number of rows of a reduced tridiagonal system after MPI_Gather
    
    nr_rd = 2;
    nr_rt = nr_rd*nprocs;

    plan.myrank = myrank;
    plan.nprocs = nprocs;
    plan.gather_rank = gather_rank;
    plan.ptdma_world = mpi_world;
    plan.n_row_rt = nr_rt;

    plan.A_rd.resize(nr_rd), plan.B_rd.resize(nr_rd), plan.C_rd.resize(nr_rd), plan.D_rd.resize(nr_rd);
    plan.A_rt.resize(nr_rt), plan.B_rt.resize(nr_rt), plan.C_rt.resize(nr_rt), plan.D_rt.resize(nr_rt);
}

void PaScaL_TDMA::PaScaL_TDMA_plan_single_destroy(ptdma_plan_single& plan) {
    // 원래는 A, B, C, D 를 deallocate를 해야하는데 cpp에서는 굳이 필요가 없음
    // 자동으로 메모리 해제 해준다.
}

void PaScaL_TDMA::PaScaL_TDMA_single_solve(ptdma_plan_single& plan, 
                                std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& D, int n_row) {

    // 1) 초기 작업 및 reduced system 만든다.
    // A_rd, ... 여기에 저장
    A[0] = A[0]/B[0];
    D[0] = D[0]/B[0];
    C[0] = C[0]/B[0];

    A[1] = A[1]/B[1];
    D[1] = D[1]/B[1];
    C[1] = C[1]/B[1];
    
    double r;
    for (int i=2; i<n_row; ++i) {
        r = 1/(B[i] - A[i]*C[i-1]);
        D[i] = r*(D[i] - A[i]*D[i-1]);
        C[i] = r*C[i];
        A[i] = -r*A[i]*A[i-1];
    }

    // Reduction step : elimination of upper diagonal elements
    for (int i=n_row-3; i>=1; --i) {
        D[i] = D[i] - C[i]*D[i+1];
        A[i] = A[i] - C[i]*A[i+1];
        C[i] = -C[i]*C[i+1];
    }

    r = 1/(1 - A[1]*C[0]);
    D[0] = r*(D[0]-C[0]*D[1]);
    A[0] = r*A[0];
    C[0] = -r*C[0]*C[1];

    // 여기서 통신해서 푼다.

    // 2) MPI_Igather 로 reduced system 을 gather_rank 로 모은다.
    // A_rd -> A_rt 로 전송
    // A_rt, ... 여기에 저장

    // Construct a reduced tridiagonal system of equations per each rank. Each process has two reduced rows.
    plan.A_rd[0] = A[0];    plan.A_rd[1] = A[n_row-1];
    plan.B_rd[0] = 1;       plan.B_rd[1] = 1;
    plan.C_rd[0] = C[0];    plan.C_rd[1] = C[n_row-1];
    plan.D_rd[0] = D[0];    plan.D_rd[1] = D[n_row-1];

    // Gather the coefficients of the reduced tridiagonal system to a defined rank, plan%gather_rank.
    std::vector<int> request(4);

    MPI_Igather(plan.A_rd.data(), 2, MPI_DOUBLE, 
                plan.A_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[0]);
    MPI_Igather(plan.B_rd.data(), 2, MPI_DOUBLE, 
                plan.B_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[1]);
    MPI_Igather(plan.C_rd.data(), 2, MPI_DOUBLE, 
                plan.C_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[2]);
    MPI_Igather(plan.D_rd.data(), 2, MPI_DOUBLE, 
                plan.D_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[3]);
    MPI_Waitall(4, request.data(), MPI_STATUS_IGNORE);

    // A_rt 에 잘 담기고 있는지 확인 해야하는데 일단 pass (귀찮음)
    // int myrank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // std::cout << "A_rd| " << "myrank: " << myrank << "|";
    // for (int i=0; i<plan.n_row_rt; ++i) {
    //     std::cout << plan.A_rd[i] << " ";
    // }

    // 3) gather_rank가 reduced system 을 tdma_single 으로 푼다.
    if (plan.myrank==plan.gather_rank) {
        // auto cx = topo.commX();
        // std::cout << "myrank_x: " << cx.myrank << "|D_rt| ";
        // for (int i=0; i<4; ++i) {
        //     std::cout << plan.D_rt[i] << " ";
        // }
        // std::cout << " \n";
        // std::cout << "myrank_x: " << cx.myrank << "|A_rt| ";
        // for (int i=0; i<4; ++i) {
        //     std::cout << plan.A_rt[i] << " ";
        // }
        // std::cout << " \n";
        // std::cout << "myrank_x: " << cx.myrank << "|B_rt| ";
        // for (int i=0; i<4; ++i) {
        //     std::cout << plan.B_rt[i] << " ";
        // }
        // std::cout << " \n";
        // std::cout << "myrank_x: " << cx.myrank << "|C_rt| ";
        // for (int i=0; i<4; ++i) {
        //     std::cout << plan.C_rt[i] << " ";
        // }
        
        tdma_single(plan.A_rt, plan.B_rt, plan.C_rt, plan.D_rt, plan.n_row_rt);
        // auto cx = topo.commX();
        // std::cout << "myrank_x: " << cx.myrank << "|D_rt|";
        // for (int i=0; i<4; ++i) {
        //     std::cout << plan.D_rt[i] << " ";
        // }
        // std::cout << "------- \n";
    };
    
    // 4) MPI_Iscatter로 각 랭크에 값을 뿌린다.
    MPI_Iscatter(plan.D_rt.data(), 2, MPI_DOUBLE,
                plan.D_rd.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[0]);
    MPI_Waitall(1, request.data(), MPI_STATUS_IGNORE);

    // auto cx = topo.commX();
    // std::cout << "myrank_x: " << cx.myrank << "|D_rd|";
    // for (int i=0; i<2; ++i) {
    //     std::cout << plan.D_rd[i] << " ";
    // }

    // 5) 이젠 그냥 알아서 푼다.

    D[0] = plan.D_rd[0];
    D[n_row-1] = plan.D_rd[1];
    for (int i=1; i<n_row-1; ++i) {
        D[i] = D[i] - A[i]*D[0] - C[i]*D[n_row-1];
    };

    // auto cx = topo.commX();
    // std::cout << "myrank_x: " << cx.myrank << "|D|: ";
    // for (int i=0; i<n_row; ++i) {
    //     std::cout << D[i]<< " ";
    // }
}

void PaScaL_TDMA::PaScaL_TDMA_single_solve_cycle(ptdma_plan_single& plan, 
                                std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& D, int n_row) {

    // The modified Thomas algorithm : elimination of lower diagonal elements.
    A[0] = A[0]/B[0];
    D[0] = D[0]/B[0];
    C[0] = C[0]/B[0];

    A[1] = A[1]/B[1];
    D[1] = D[1]/B[1];
    C[1] = C[1]/B[1];

    double r;
    for (int i=2; i<n_row; ++i) {
        r = 1/(B[i] - A[i]*C[i-1]);
        D[i] = r*(D[i] - A[i]*D[i-1]);
        C[i] = r*C[i];
        A[i] = -r*A[i]*A[i-1];
    }

    // The modified Thomas algorithm : elimination of upper diagonal elements.
    for (int i=n_row-3; i>=1; --i) {
        D[i] = D[i] - C[i]*D[i+1];
        A[i] = A[i] - C[i]*A[i+1];
        C[i] = -C[i]*C[i+1];
    }

    r = 1/(1 - A[1]*C[0]);
    D[0] = r*(D[0]-C[0]*D[1]);
    A[0] = r*A[0];
    C[0] = -r*C[0]*C[1];

    // Construct a reduced tridiagonal system of equations per each rank. Each process has two reduced rows.
    plan.A_rd[0] = A[0];    plan.A_rd[1] = A[n_row-1];
    plan.B_rd[0] = 1;       plan.B_rd[1] = 1;
    plan.C_rd[0] = C[0];    plan.C_rd[1] = C[n_row-1];
    plan.D_rd[0] = D[0];    plan.D_rd[1] = D[n_row-1];

    // Gather the coefficients of the reduced tridiagonal system to a defined rank, plan%gather_rank.
    std::vector<int> request(4);

    MPI_Igather(plan.A_rd.data(), 2, MPI_DOUBLE, 
                plan.A_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[0]);
    MPI_Igather(plan.B_rd.data(), 2, MPI_DOUBLE, 
                plan.B_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[1]);
    MPI_Igather(plan.C_rd.data(), 2, MPI_DOUBLE, 
                plan.C_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[2]);
    MPI_Igather(plan.D_rd.data(), 2, MPI_DOUBLE, 
                plan.D_rt.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[3]);
    MPI_Waitall(4, request.data(), MPI_STATUS_IGNORE);

    // Solve the reduced cyclic tridiagonal system on plan%gather_rank.
    if (plan.myrank==plan.gather_rank) {
        tdma_cycl_single(plan.A_rt, plan.B_rt, plan.C_rt, plan.D_rt, plan.n_row_rt);
    };

    // Distribute the solutions to each rank.
    MPI_Iscatter(plan.D_rt.data(), 2, MPI_DOUBLE,
                plan.D_rd.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[0]);
                
    MPI_Waitall(1, request.data(), MPI_STATUS_IGNORE);

    // Update solutions of the modified tridiagonal system with the solutions of the reduced tridiagonal system.
    D[0] = plan.D_rd[0];
    D[n_row-1] = plan.D_rd[1];
    for (int i=1; i<n_row-1; ++i) {
        D[i] = D[i] - A[i]*D[0] - C[i]*D[n_row-1];
    };
}


