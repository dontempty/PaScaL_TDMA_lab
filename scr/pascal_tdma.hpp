#include <mpi.h>
#include <vector>

class PaScaL_TDMA {
public:
    struct ptdma_plan_single {
        MPI_Comm ptdma_world;
        int n_row_rt;
        int gather_rank, myrank, nprocs;
        std::vector<double> A_rd, B_rd, C_rd, D_rd;
        std::vector<double> A_rt, B_rt, C_rt, D_rt;
    };

    void PaScaL_TDMA_plan_single_create(ptdma_plan_single& plan,
                                        int myrank, int nprocs,
                                        MPI_Comm mpi_world,
                                        int gather_rank);
    void PaScaL_TDMA_plan_single_destroy(ptdma_plan_single& plan);
};
