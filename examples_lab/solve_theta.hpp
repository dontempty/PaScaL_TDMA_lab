// solve_theta.hpp
#ifndef SOLVE_THETA_HPP
#define SOLVE_THETA_HPP

#include <mpi.h>
#include <vector>
#include "mpi_topology.hpp"
#include "mpi_subdomain.hpp"
#include "../src/pascal_tdma.hpp"

class solve_theta {
public:
    solve_theta(const GlobalParams& params,
                const MPITopology& topo,
                MPISubdomain& sub);

    // rankx, npx 등 필요한 파라미터를 인자로 받아야 합니다
    void solve_theta_plan_single(std::vector<double>& theta);

private:
    const GlobalParams& params;
    const MPITopology& topo;
    MPISubdomain& sub;                         
};

#endif // SOLVE_THETA_HPP