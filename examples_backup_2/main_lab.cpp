// #include "headers.hpp"
#include <mpi.h>
#include "global.hpp"
#include "mpi_topology.hpp"
#include "mpi_subdomain.hpp"

#include <iostream>
#include <array>
#include <string>

int main(int argc, char** argv) {

    int nprocs, myrank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    // 1) 전역 파라미터 읽기
    GlobalParams params;
    params.load(argv[1]); 
    // std::cout << "Prandtl number: " << params.Pr << "\n";
    // std::cout << "Grid size: nx=" << params.nx
    //           << ", ny=" << params.ny
    //           << ", nz=" << params.nz << "\n";
    // std::cout << "MPI process refine: ("
    //           << params.np_dim[0] << ", "
    //           << params.np_dim[1] << ", "
    //           << params.np_dim[2] << ")\n";
    
    // 2) Cartesian 토폴로지 생성
    MPITopology topo({params.np_dim[0], params.np_dim[1], params.np_dim[2]}, {true, false, true});
    topo.make();
    auto cx = topo.commX();
    auto cy = topo.commY();
    auto cz = topo.commZ();
    // std::cout <<"Global 3D-cart rank: " << myrank
    //           << " / X-subcomm rank: " << cx.myrank
    //           << " / nprocs: "      << cx.nprocs
    //           << " / west: "      << cx.west_rank
    //           << " / east: "      << cx.east_rank
    //           << "\n";

    
    // 3) npx, rankx 등 변수 정의
    int npx   = params.np_dim[0], rankx = cx.myrank;
    int npy   = params.np_dim[1], ranky = cy.myrank;
    int npz   = params.np_dim[2], rankz = cz.myrank;

    // 4) 서브도메인 객체 만들기
    MPISubdomain sub;
    sub.make(params, npx, rankx, npy, ranky, npz, rankz);
    // // ── 디버그1: 서브도메인 크기/인덱스 확인
    std::cout << "Rank XYZ = " << rankx << ranky << rankz << "|"
              << "nx_sub=" << sub.nx_sub
              << ", ny_sub=" << sub.ny_sub
              << ", nz_sub=" << sub.nz_sub << "\n    "
              << ", ista..iend = "
              << sub.ista << " .. " << sub.iend
              << ", jsta..jend = "
              << sub.jsta << " .. " << sub.jend
              << ", ksta..kend = "
              << sub.ksta << " .. " << sub.kend
              << "\n" << std::endl;

    MPI_Barrier(MPI_COMM_WORLD);

    // 5) ghostcell 통신용 MPI_Datatype 생성
    sub.makeGhostcellDDType();

    // 6) 메쉬(좌표, dx·dy·dz) 계산
    sub.mesh(params, rankx, ranky, rankz, npx, npy, npz);
    // // ── 디버그2: 좌표 벡터 일부 출력
    // std::cout << "[Debug] x_sub[0..] = ";
    // for (int i = 0; i <= std::min(sub.nx_sub, 5); ++i)
    //     std::cout << sub.x_sub[i] << " ";
    // std::cout << "...\n";
    // std::cout << "[Debug] y_sub[0..] = ";
    // for (int j = 0; j <= std::min(sub.ny_sub, 5); ++j)
    //     std::cout << sub.y_sub[j] << " ";
    // std::cout << "...\n";
    // std::cout << "[Debug] z_sub[0..] = ";
    // for (int k = 0; k <= std::min(sub.nz_sub, 5); ++k)
    //     std::cout << sub.z_sub[k] << " ";
    // std::cout << "...\n";

    // 7) y-방향 경계 인덱스 설정
    sub.indices(params, ranky, npy);

    // 8) theta 필드 (subdomain 크기에 맞춘 flat 배열) 준비
    int nx1 = sub.nx_sub + 1;
    int ny1 = sub.ny_sub + 1;
    int nz1 = sub.nz_sub + 1;
    std::vector<double> theta(nx1 * ny1 * nz1, 0.0);

    // 9) 내부 값 초기화
    sub.initialization(theta.data(), params, ranky, npy, myrank);

    // 10) 경계 조건 추출
    sub.boundary(theta.data(), params, ranky, npy);
    // // ── 디버그3: 경계 버퍼 일부 확인
    // std::cout << "[Debug] thetaBC3_sub[0..] = ";
    // for (int idx = 0; idx < std::min((int)sub.thetaBC3_sub.size(), 6); ++idx)
    //     std::cout << sub.thetaBC3_sub[idx] << " ";
    // std::cout << "...\n";

    // 11) ghostcell 교환
    sub.ghostcellUpdate(theta.data(), cx, cy, cz, params);
    std::cout << "myrank: " << myrank << "|"
              << "Rank XYZ = " << rankx << ranky << rankz << "|"
              << "n_sub = " << nx1 << ny1 << nz1 << "|";
    for (int i=0; i<nx1 * ny1 * nz1; ++i) {
        std::cout << theta[i] << " ";
    }
    std::cout << "\n";

    // 12) 정리
    sub.clean();
    topo.clean();
    MPI_Finalize();

    return 0;
}
