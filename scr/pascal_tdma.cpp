#include <mpi.h>
#include <vector>
#include "pascal_tdma.hpp"

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
                                std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d) {
    
    // 여기서 통신해서 푼다.

    // 1) 초기 작업 및 reduced system 만든다.
    // A_rd, ... 여기에 저장

    // 2) MPI_Igather 로 reduced system 을 gather_rank 로 모은다.
    // A_rd -> A_rt 로 전송
    // A_rt, ... 여기에 저장

    // 3) gather_rank가 reduced system 을 tdma_single 으로 푼다.
    // if(plan%myrank == plan.gather_rank) then

    // 4) MPI_Iscatter로 각 랭크에 값을 뿌린다.

    // 5) 이젠 그냥 알아서 푼다.
}