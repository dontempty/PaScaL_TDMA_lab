#include <mpi.h>
#include <vector>
#include "iostream"
#include <numeric>

#include "pascal_tdma.hpp"
#include "tdmas.hpp"
#include "../examples_lab/mpi_subdomain.hpp"
#include "../examples_lab/mpi_topology.hpp"
#include "para_range.hpp"

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

    // 3) gather_rank가 reduced system 을 tdma_single 으로 푼다.
    if (plan.myrank==plan.gather_rank) {
        tdma_single(plan.A_rt, plan.B_rt, plan.C_rt, plan.D_rt, plan.n_row_rt);
    };
    
    // 4) MPI_Iscatter로 각 랭크에 값을 뿌린다.
    MPI_Iscatter(plan.D_rt.data(), 2, MPI_DOUBLE,
                plan.D_rd.data(), 2, MPI_DOUBLE,
                plan.gather_rank, plan.ptdma_world, &request[0]);
    MPI_Waitall(1, request.data(), MPI_STATUS_IGNORE);

    // 5) 이젠 그냥 알아서 푼다.
    D[0] = plan.D_rd[0];
    D[n_row-1] = plan.D_rd[1];
    for (int i=1; i<n_row-1; ++i) {
        D[i] = D[i] - A[i]*D[0] - C[i]*D[n_row-1];
    };
}

void PaScaL_TDMA::PaScaL_TDMA_single_solve_cycle(ptdma_plan_single& plan, 
                                std::vector<double>& A, std::vector<double>& B, std::vector<double>& C, std::vector<double>& D,
                                int n_row) {

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

// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
// --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

// Create a plan for a single tridiagonal system of equations.
void PaScaL_TDMA::PaScaL_TDMA_plan_many_create(ptdma_plan_many& plan, int n_sys, int myrank, int nprocs, MPI_Comm mpi_world) {

    int i, ierr;
    int ista, iend;                                         // First and last indices of assigned range in many tridiagonal systems of equations 
    std::vector<int> bigsize(2), subsize(2), start(2);      // Temporary variables of derived data type (DDT)
    int ns_rd, nr_rd;                                       // Dimensions of many reduced tridiagonal systems
    int ns_rt, nr_rt;                                       // Dimensions of many reduced tridiagonal systems after transpose
    std::vector<int> ns_rt_array(nprocs);                   // Array specifying the number of tridiagonal systems for each process after transpose

    // [ ][ ] | [ ][ ] | [ ][ ]
    // [ ][ ] | [ ][ ] | [ ][ ]
    // [ ][ ] | [ ][ ] | [ ][ ]
    // ------------------------
    // [ ][ ] | [ ][ ] | [ ][ ]
    // [ ][ ] | [ ][ ] | [ ][ ]
    // [0][0] | [0][0] | [0][0]

    // n_sys : x 축으로 푼다고 가정할 때, 한 불록의 y 개수 (3)
    // 각 y 마다 reduced system을 가지는데, 한 블록 기준으로 y에 있는거 까지 다 모은다.

    plan.nprocs = nprocs;

    // Specify dimensions for reduced systems.
    ns_rd = n_sys;
    nr_rd = 2;

    // Specify dimensions for reduced systems after transpose.
    // ns_rt         : divide the number of tridiagonal systems of equations per each process  
    // ns_rt_array   : save the ns_rt in ns_rt_array for defining the DDT
    // nr_rt         : dimensions of the reduced tridiagonal systems in the solving direction, nr_rd*nprocs
    para_range(1, ns_rd, nprocs, myrank, ista, iend);
    ns_rt = iend - ista + 1;
    MPI_Allgather(&ns_rt, 1, MPI_INT,
                  ns_rt_array.data(), 1, MPI_INT,
                  mpi_world);
    nr_rt = nr_rd*nprocs;

    // Assign plan variables and allocate coefficient arrays.
    plan.n_sys_rt = ns_rt;
    plan.n_row_rt = nr_rt;
    plan.ptdma_world = mpi_world;

    // 2d array (ns_rd, nr_rd) or (ns_rt, nr_rt)
    // 여기서는 init이 없기 때문인지 resize를 해주는게 더 빠름. (re declare랑 비교했을 때)
    plan.A_rd.resize(ns_rd);
    plan.B_rd.resize(ns_rd);
    plan.C_rd.resize(ns_rd);
    plan.D_rd.resize(ns_rd);
    for (int i = 0; i < ns_rd; ++i) {
        plan.A_rd[i].resize(nr_rd);
        plan.B_rd[i].resize(nr_rd);
        plan.C_rd[i].resize(nr_rd);
        plan.D_rd[i].resize(nr_rd);
    }

    plan.A_rt.resize(ns_rt);
    plan.B_rt.resize(ns_rt);
    plan.C_rt.resize(ns_rt);
    plan.D_rt.resize(ns_rt);
    for (int i = 0; i < ns_rt; ++i) {
        plan.A_rt[i].resize(nr_rt);
        plan.B_rt[i].resize(nr_rt);
        plan.C_rt[i].resize(nr_rt);
        plan.D_rt[i].resize(nr_rt);
    }

    // Building the DDTs.
    plan.ddtype_Fs.resize(nprocs), plan.ddtype_Bs.resize(nprocs);
    int sum = 0;
    for (i=0; i<nprocs; ++i) {
        // DDT for sending coefficients of the reduced tridiagonal systems using MPI_Ialltoallw communication.
        bigsize[0] = ns_rd;
        bigsize[1] = nr_rd;
        subsize[0] = ns_rt_array[i];
        subsize[1] = nr_rd;
        start[0] = sum;
        start[1] = 0;
        MPI_Type_create_subarray(2, bigsize.data(), subsize.data(), start.data(),
                                        MPI_ORDER_C, MPI_DOUBLE,
                                        &plan.ddtype_Fs[i]);
        MPI_Type_commit(&plan.ddtype_Fs[i]);
        
        // DDT for receiving coefficients for the transposed systems of reduction using MPI_Ialltoallw communication.
        bigsize[0] = ns_rt;
        bigsize[1] = nr_rt;
        subsize[0] = ns_rt;
        subsize[1] = nr_rd;
        start[0] = 0;
        start[1] = nr_rd*i;
        MPI_Type_create_subarray(2, bigsize.data(), subsize.data(), start.data(), 
                                        MPI_ORDER_C, MPI_DOUBLE,
                                        &plan.ddtype_Bs[i]);
        MPI_Type_commit(&plan.ddtype_Bs[i]);
    }

    // Buffer counts and displacements for MPI_Ialltoallw.
    // All buffer counts are 1 and displacements are 0 due to the defined DDT.
    // resize + init 이랑 re declare 성능을 비교 했는데 후자가 조금 더 빠름 (init도 같이 진행하기 때문)
    plan.count_send = std::vector<int>(nprocs, 1);
    plan.displ_send = std::vector<int>(nprocs, 0);
    plan.count_recv = std::vector<int>(nprocs, 1);
    plan.displ_recv = std::vector<int>(nprocs, 0);     

    // Deallocate local array.
    if (!ns_rt_array.empty()) {
        ns_rt_array.clear();         // 내용 제거
        // ns_rt_array.shrink_to_fit(); // capacity까지 제거
    }
}

void PaScaL_TDMA::PaScaL_TDMA_plan_many_destroy(ptdma_plan_many& plan, int nprocs) {

    for (int i=0; i<nprocs; ++i) {
        MPI_Type_free(&plan.ddtype_Fs[i]);
        MPI_Type_free(&plan.ddtype_Bs[i]);
    }

    plan.ddtype_Fs.shrink_to_fit(); plan.ddtype_Fs.shrink_to_fit();
    plan. count_send.shrink_to_fit(); plan.displ_send.shrink_to_fit();
    plan.count_recv.shrink_to_fit(); plan.displ_recv.shrink_to_fit();
    plan.A_rd.shrink_to_fit(); plan.B_rd.shrink_to_fit(); plan.C_rd.shrink_to_fit(); plan.D_rd.shrink_to_fit();
    plan.A_rt.shrink_to_fit(); plan.B_rt.shrink_to_fit(); plan.C_rt.shrink_to_fit(); plan.D_rt.shrink_to_fit();
}

void PaScaL_TDMA::PaScaL_TDMA_many_solve(ptdma_plan_many& plan,
                                std::vector<std::vector<double>>& A, 
                                std::vector<std::vector<double>>& B, 
                                std::vector<std::vector<double>>& C, 
                                std::vector<std::vector<double>>& D,
                                int n_sys, int n_row) {

    // Temporary variables for computation and parameters for MPI functions.
    int i, j;
    std::vector<MPI_Request> request(4);
    double r;
    int idx;

    for (i = 0; i < n_sys; ++i) {
        A[i][0] /= B[i][0];
        D[i][0] /= B[i][0];
        C[i][0] /= B[i][0];

        A[i][1] /= B[i][1];
        D[i][1] /= B[i][1];
        C[i][1] /= B[i][1];
    }

    for (j = 2; j < n_row; ++j) {
        for (i = 0; i < n_sys; ++i) {
            r = 1.0 / (B[i][j] - A[i][j] * C[i][j - 1]);
            D[i][j] = r * (D[i][j] - A[i][j] * D[i][j - 1]);
            C[i][j] = r * C[i][j];
            A[i][j] = -r * A[i][j] * A[i][j - 1];
        }
    }

    for (j = n_row - 3; j >= 1; --j) {
        for (i = 0; i < n_sys; ++i) {
            D[i][j] -= C[i][j] * D[i][j + 1];
            A[i][j] -= C[i][j] * A[i][j + 1];
            C[i][j] = -C[i][j] * C[i][j + 1];
        }
    }

    for (i = 0; i < n_sys; ++i) {
        r = 1.0 / (1.0 - A[i][1] * C[i][0]);
        D[i][0] = r * (D[i][0] - C[i][0] * D[i][1]);
        A[i][0] *= r;
        C[i][0] = -r * C[i][0] * C[i][1];

        plan.A_rd[i][0] = A[i][0];
        plan.A_rd[i][1] = A[i][n_row - 1];
        plan.B_rd[i][0] = 1.0;
        plan.B_rd[i][1] = 1.0;
        plan.C_rd[i][0] = C[i][0];
        plan.C_rd[i][1] = C[i][n_row - 1];
        plan.D_rd[i][0] = D[i][0];
        plan.D_rd[i][1] = D[i][n_row - 1];
    }

    std::vector<double> flat_A_rd(n_sys * 2);
    std::vector<double> flat_B_rd(n_sys * 2);
    std::vector<double> flat_C_rd(n_sys * 2);
    std::vector<double> flat_D_rd(n_sys * 2);
    std::vector<double> flat_A_rt(plan.n_sys_rt * plan.n_row_rt);
    std::vector<double> flat_B_rt(plan.n_sys_rt * plan.n_row_rt);
    std::vector<double> flat_C_rt(plan.n_sys_rt * plan.n_row_rt);
    std::vector<double> flat_D_rt(plan.n_sys_rt * plan.n_row_rt);

    for (int i = 0; i < n_sys; ++i) {
        flat_A_rd[i * 2 + 0] = plan.A_rd[i][0];
        flat_A_rd[i * 2 + 1] = plan.A_rd[i][1];

        flat_B_rd[i * 2 + 0] = plan.B_rd[i][0];
        flat_B_rd[i * 2 + 1] = plan.B_rd[i][1];

        flat_C_rd[i * 2 + 0] = plan.C_rd[i][0];
        flat_C_rd[i * 2 + 1] = plan.C_rd[i][1];

        flat_D_rd[i * 2 + 0] = plan.D_rd[i][0];
        flat_D_rd[i * 2 + 1] = plan.D_rd[i][1];
    }
    
    MPI_Ialltoallw(flat_A_rd.data(), plan.count_send.data(), plan.displ_send.data(), plan.ddtype_Fs.data(),
                    flat_A_rt.data(), plan.count_recv.data(), plan.displ_recv.data(), plan.ddtype_Bs.data(),
                    plan.ptdma_world, &request[0]);
    MPI_Ialltoallw(flat_B_rd.data(), plan.count_send.data(), plan.displ_send.data(), plan.ddtype_Fs.data(),
                    flat_B_rt.data(), plan.count_recv.data(), plan.displ_recv.data(), plan.ddtype_Bs.data(),
                    plan.ptdma_world, &request[1]);
    MPI_Ialltoallw(flat_C_rd.data(), plan.count_send.data(), plan.displ_send.data(), plan.ddtype_Fs.data(),
                    flat_C_rt.data(), plan.count_recv.data(), plan.displ_recv.data(), plan.ddtype_Bs.data(),
                    plan.ptdma_world, &request[2]);     
    MPI_Ialltoallw(flat_D_rd.data(), plan.count_send.data(), plan.displ_send.data(), plan.ddtype_Fs.data(),
                    flat_D_rt.data(), plan.count_recv.data(), plan.displ_recv.data(), plan.ddtype_Bs.data(),
                    plan.ptdma_world, &request[3]);
   
    MPI_Waitall(4, request.data(), MPI_STATUSES_IGNORE);

    // 각 랭크에서 통신 직후 확인용 출력
    // int myrank;
    // MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    // std::cout << "=== Rank " << myrank << " | flat_D_rt ===\n";
    // for (int i = 0; i < plan.n_sys_rt; ++i) {
    //     std::cout << "  sys " << i << ": \n";
    //     for (int j = 0; j < plan.n_row_rt; ++j) {
    //         std::cout << flat_D_rt[i * plan.n_row_rt + j] << " ";
    //     }
    //     std::cout << "\n";
    // }

    for (int i=0; i<plan.n_sys_rt; ++i) {
        for (int j=0; j<plan.n_row_rt; ++j) {
            plan.A_rt[i][j] = flat_A_rt[i * plan.n_row_rt + j];
            plan.B_rt[i][j] = flat_B_rt[i * plan.n_row_rt + j];
            plan.C_rt[i][j] = flat_C_rt[i * plan.n_row_rt + j];
            plan.D_rt[i][j] = flat_D_rt[i * plan.n_row_rt + j];
        }
    }

    tdma_many(plan.A_rt, plan.B_rt, plan.C_rt, plan.D_rt, plan.n_sys_rt, plan.n_row_rt);

    for (int i=0; i<plan.n_sys_rt; ++i) {
        for (int j=0; j<plan.n_row_rt; ++j) {
            flat_D_rt[i * plan.n_row_rt + j] = plan.D_rt[i][j];
        }
    }

    // 일단 이렇게 하면 작동은 잘 되는데 왜 잘 작동하는지는 모르겠다
    // D_rt -> D_rd의 경우 보내는 입장이 바뀌었기 때문에 ddtype을 (ddtype_Bs, ddtype_F) 이렇게 사용해야 하는줄 알았는데 반대로 사용해야 잘 작동함.
    MPI_Ialltoallw(flat_D_rt.data(), plan.count_recv.data(), plan.displ_recv.data(), plan.ddtype_Fs.data(),
               flat_D_rd.data(), plan.count_send.data(), plan.displ_send.data(), plan.ddtype_Bs.data(),
               plan.ptdma_world, &request[0]);

    MPI_Waitall(1, request.data(), MPI_STATUSES_IGNORE);

    for (int i = 0; i < n_sys; ++i) {
        plan.D_rd[i][0] = flat_D_rd[i * 2 + 0];
        plan.D_rd[i][1] = flat_D_rd[i * 2 + 1];
    }

    for (i = 0; i < n_sys; ++i) {
        D[i][0] = plan.D_rd[i][0];
        D[i][n_row - 1] = plan.D_rd[i][1];
    }

    for (j = 1; j < n_row - 1; ++j) {
        for (i = 0; i < n_sys; ++i) {
            D[i][j] -= A[i][j] * D[i][0] + C[i][j] * D[i][n_row - 1];
        }
    }
}