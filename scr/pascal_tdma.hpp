#include <mpi.h>
#include <vector>

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

    void PaScaL_TDMA_plan_single_create(ptdma_plan_single& plan,
                                        int myrank, int nprocs,
                                        MPI_Comm mpi_world,
                                        int gather_rank);
                                        
    void PaScaL_TDMA_plan_single_destroy(ptdma_plan_single& plan);

    void PaScaL_TDMA_single_solve(ptdma_plan_single& plan, 
                                  std::vector<double>& a, std::vector<double>& b, std::vector<double>& c, std::vector<double>& d);
};
