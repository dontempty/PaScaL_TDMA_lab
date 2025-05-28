// mpi_topology.hpp
#ifndef MPI_TOPOLOGY_HPP
#define MPI_TOPOLOGY_HPP

#include <mpi.h>
#include <array>
#include <stdexcept>

struct CartComm1D {
    int      myrank, nprocs, west_rank, east_rank;
    MPI_Comm comm;
};

class MPITopology {
public:
    MPITopology(const std::array<int,2>& dims,
                const std::array<bool,2>& periods);
    void make();
    void clean();
    MPI_Comm            worldCart() const;
    const CartComm1D&   commX()     const;
    const CartComm1D&   commY()     const;
    // const CartComm1D&   commZ()     const;
private:
    void defineSubcomm(int dim, CartComm1D& sub);
    std::array<int,2>   dims_;
    std::array<bool,2>  periods_;
    MPI_Comm            world_cart_;
    CartComm1D          comm_x_, comm_y_, comm_z_;
};

#endif // MPI_TOPOLOGY_HPP
