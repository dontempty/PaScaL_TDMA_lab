// #include "headers.hpp"
#include <mpi.h>
#include "global.hpp"         // params 생성
#include "mpi_topology.hpp"   // topo 생성
#include "mpi_subdomain.hpp"  // MPISubdomain sub; 그냥 extern 함
#include "solve_theta.hpp"

#include <iostream>
#include <array>
#include <string>

int main(int argc, char** argv) {

    int nprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    params.load(argv[1]); 

    // x, y, z
    topo.init(
      { params.np_dim[0], params.np_dim[1] },
      { false, false }
    );

    topo.make();

    auto cx = topo.commX();
    auto cy = topo.commY();

    // 3) npx, rankx 등 변수 정의
    int npx   = params.np_dim[0], rankx = cx.myrank;
    int npy   = params.np_dim[1], ranky = cy.myrank;

    // 4) 서브도메인 객체 만들기
    sub.make(params, npx, rankx, npy, ranky);

    // 5) ghostcell 통신용 MPI_Datatype 생성
    sub.makeGhostcellDDType();

    sub.mesh(params, rankx, ranky, npx, npy);

    // 7) x, y-방향 경계 인덱스 설정
    sub.indices(params, rankx, npx, ranky, npy);

    // 8) theta 필드 (subdomain 크기에 맞춘 flat 배열) 준비
    std::vector<double> theta((sub.ny_sub + 1) * (sub.nx_sub + 1), 0.0);

    // 9) 내부 값 초기화
    // initialization(theta.data(), params, ranky, npy);
    sub.initialization(theta.data(), params);
    MPI_Barrier(MPI_COMM_WORLD);

    // 10) ghostcell 교환
    sub.ghostcellUpdate(theta.data(), cx, cy, params);

    // 11) bdy값 저장하기
    sub.boundary(theta.data(), params, rankx, npx, ranky, npy);

    // 12) 솔버 호출해서 실행하기
    solve_theta solver;
    solver.solve_theta_plan_single(theta);
    
    sub.clean();
    topo.clean();
    MPI_Finalize();
};
