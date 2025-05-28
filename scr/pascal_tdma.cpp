#include <mpi.h>
#include <vector>
#include "pascal_tdma.hpp"

class PaScaL_TDMA {

public:
    struct ptdma_plan_single {
        int ptdma_world;    // Single dimensional subcommunicator to assemble data for the reduced TDMA
        int n_row_rt;       // Number of rows of a reduced tridiagonal system after MPI_Gather
        
        int gather_rank;    // Destination rank of MPI_Igather
        int myrank;         // Current rank ID in the communicator of ptdma_world
        int nprocs;         // Communicator size of ptdma_world
        
        // Coefficient arrays after reduction, a: lower, b: diagonal, c: upper, d: rhs.
        std::vector<double> A_rd, B_rd, C_rd, D_rd;

        // Coefficient arrays after transpose of a reduced system, a: lower, b: diagonal, c: upper, d: rhs
        std::vector<double> A_rt, B_rt, C_rt, D_rt;

    };

    // Create a plan for a single tridiagonal system of equations.
    void PaScaL_TDMA_plan_single_create(ptdma_plan_single& plan, int myrank, int nprocs, MPI_Comm mpi_world, int gather_rank) {

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

    void PaScaL_TDMA_plan_single_destroy(ptdma_plan_single& plan) {
        // 원래는 A, B, C, D 를 deallocate를 해야하는데 cpp에서는 굳이 필요가 없음
        // 자동으로 메모리 해제 해준다.
    }
};