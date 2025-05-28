// global_params.hpp
#ifndef GLOBAL_PARAMS_HPP
#define GLOBAL_PARAMS_HPP

#include <array>
#include <string>

#ifdef USE_MPI
#include <mpi.h>
#endif

#ifdef _CUDA
#define THREAD_PARAMS
#endif

class GlobalParams {
public:
    GlobalParams() = default;

    // input file로부터 파라미터 읽기
    void load(const std::string& filename);

    // 물리 및 수치 파라미터
    double Pr, Ra;
    int    Tmax;

    int nx, ny, nz;
    int nxm, nym, nzm;
    int nxp, nyp, nzp;

    double dt, dtStart, tStart;
    double lx, ly, lz;
    double dx, dy, dz;

    double theta_cold, theta_hot, alphaG, nu, Ct;

#ifdef THREAD_PARAMS
    // CUDA 쓰레드 관련 파라미터
    int thread_in_x, thread_in_y, thread_in_z;
    int thread_in_x_pascal, thread_in_y_pascal;
#endif

    // MPI 프로세스 분할
    std::array<int,3> np_dim;
};

#endif // GLOBAL_PARAMS_HPP
