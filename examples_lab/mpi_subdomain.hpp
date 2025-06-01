// mpi_subdomain.hpp
#ifndef MPI_SUBDOMAIN_HPP
#define MPI_SUBDOMAIN_HPP

#include <mpi.h>
#include <vector>
#include <array>
#include "global.hpp"
#include "mpi_topology.hpp"

extern const double Pi;

class MPISubdomain {
public:
    // 기본 생성자: MPI_Datatype을 MPI_DATATYPE_NULL로 초기화
    MPISubdomain()
      : ddtype_sendto_E(MPI_DATATYPE_NULL), ddtype_recvfrom_W(MPI_DATATYPE_NULL),
        ddtype_sendto_W(MPI_DATATYPE_NULL), ddtype_recvfrom_E(MPI_DATATYPE_NULL),
        ddtype_sendto_N(MPI_DATATYPE_NULL), ddtype_recvfrom_S(MPI_DATATYPE_NULL),
        ddtype_sendto_S(MPI_DATATYPE_NULL), ddtype_recvfrom_N(MPI_DATATYPE_NULL) {}

    // Subdomain setup: pass simulation params and ranks
    void make(const GlobalParams& params,
              int npx, int rankx,
              int npy, int ranky);
    void clean();

    // Ghost-cell MPI derived types
    void makeGhostcellDDType();
    // Exchange ghost cells for a flat theta array
    void ghostcellUpdate(double* theta,
                         const CartComm1D& cx,
                         const CartComm1D& cy,
                         const GlobalParams& params);

    // Boundary indices in y-direction
    void indices(const GlobalParams& params,
                 int rankx, int npx,
                 int ranky, int npy);

    // Mesh coordinates and spacings
    void mesh(const GlobalParams& params,
              int rankx, int ranky,
              int npx, int npy);

    // Initialize field in subdomain
    void initialization(double* theta,
                        const GlobalParams& params);

    // Initialize field in subdomain -- debug
    void initialization_debug(double* theta,
                                  const GlobalParams& params,
                                  int myrank);

    // Boundary extraction
    void boundary(const double* theta,
                  const GlobalParams& params,
                  int rankx, int npx,
                  int ranky, int npy);

    // Subdomain sizes and offsets
    int nx_sub, ny_sub;
    int ista, iend, jsta, jend;

    // Coordinates and spacings
    std::vector<double> x_sub, y_sub;
    std::vector<double> dmx_sub, dmy_sub;

    // Ghost-cell boundary buffers (flattened)
    std::vector<double> theta_x_left_sub, theta_x_right_sub;
    std::vector<double> theta_y_left_sub, theta_y_right_sub;
    std::vector<int>    theta_x_left_index, theta_x_right_index;
    std::vector<int>    theta_y_left_index, theta_y_right_index;

    // MPI derived datatypes
    MPI_Datatype ddtype_sendto_E, ddtype_recvfrom_W, ddtype_sendto_W, ddtype_recvfrom_E;
    MPI_Datatype ddtype_sendto_N, ddtype_recvfrom_S, ddtype_sendto_S, ddtype_recvfrom_N;
};

extern MPISubdomain sub;

#endif // MPI_SUBDOMAIN_HPP