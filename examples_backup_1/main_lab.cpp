#include "headers.hpp"
#include <mpi.h>
#include "mpi_topology.hpp"

#include <iostream>
#include <array>
#include <string>

int main(int argc, char** argv) {

    int nprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    std::array<int, 3> np_dim;

    try {
        // std::cout << "[Main] The main simulation starts!" << std::endl;
        global_inputpara(argv[1], np_dim);
    } catch (const std::exception& e) {
        std::cerr << "[Exception] " << e.what() << std::endl;
        return 1;
    }
    // std::cout << "[Rank " << myrank << "] Hello from rank!" << std::endl;
    // std::cout << "===== Debug Global Input =====\n";
    // std::cout << "nx = " << nx << ", ny = " << ny << ", nz = " << nz << '\n';
    // std::cout << "Tmax = " << Tmax << '\n';
    // std::cout << "np_dim = (" << np_dim[0] << ", " << np_dim[1] << ", " << np_dim[2] << ")\n";
    // std::cout << "Pr = " << Pr << ", Ra = " << Ra << ", nu = " << nu << ", Ct = " << Ct << '\n';
    // std::cout << "==============================" << std::endl;
    
    MPITopology topo({np_dim[0], np_dim[1], np_dim[2]}, {true, false, true});
    topo.make();

    auto cx = topo.commX();
    auto cy = topo.commY();
    auto cz = topo.commZ();

    std::cout <<"Global 3D-cart rank: " << myrank
              << " / X-subcomm rank: " << cx.myrank
              << " / size: "      << cx.nprocs
              << " / west: "      << cx.west_rank
              << " / east: "      << cx.east_rank
              << "\n";

    topo.clean();
    MPI_Finalize();

    return 0;
}
