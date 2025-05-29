// mpi_subdomain.cpp
#include "mpi_subdomain.hpp"
#include <cmath>
#include "iostream"
#include <algorithm>

static void para_range(int start, int end, int nproc, int rank, int &sta, int &iend) {
    // 전체 도메인에서 ghost cell 포함한 index를 기준으로 
    // ghost cell을 뺀 인덱스의 범위 (start <= idx <= end)
    // ex: |--*--|--*--|--*--|--*--| = 실제 도메인
    // (--0--)|--start--|--2--|--3--|--end--|(--5--)  = ghost cell을 추가한 도메인
    int len = end - start + 1;
    int base = len / nproc;
    int rem  = len % nproc;
    sta = start + rank * base + std::min(rank, rem);
    iend = sta + base - 1 + (rank < rem ? 1 : 0);
}

void MPISubdomain::make(const GlobalParams& params,
                        int npx, int rankx,
                        int npy, int ranky) {

    // (--0--)|--ista--|--2--|--3--|--iend--|(--nx_sub--) 이게 각 block이 가지는 도메인이라고 가정
    para_range(1, params.nx-1, npx, rankx, ista, iend);
    nx_sub = iend - ista + 2;
    para_range(1, params.ny-1, npy, ranky, jsta, jend);
    ny_sub = jend - jsta + 2;

    // ghost cell 포함한 cell을 기준으로 하는 저장공간
    x_sub.resize(nx_sub+1);
    dmx_sub.resize(nx_sub+1);
    y_sub.resize(ny_sub+1);
    dmy_sub.resize(ny_sub+1);

    theta_x_left_sub.assign((ny_sub+1), 0.0);
    theta_x_right_sub.assign((ny_sub+1), 0.0);
    theta_y_left_sub.assign((nx_sub+1), 0.0);
    theta_y_right_sub.assign((nx_sub+1), 0.0);

    theta_x_left_index.assign(ny_sub+1, 1);
    theta_x_right_index.assign(ny_sub+1, 1);
    theta_y_left_index.assign(nx_sub+1, 1);
    theta_y_right_index.assign(nx_sub+1, 1);
}

void MPISubdomain::clean() {
    auto freeType = [](MPI_Datatype& dt) {
        if (dt != MPI_DATATYPE_NULL) MPI_Type_free(&dt);
    };
    freeType(ddtype_sendto_E);
    freeType(ddtype_recvfrom_W);
    freeType(ddtype_sendto_W);
    freeType(ddtype_recvfrom_E);
    freeType(ddtype_sendto_N);
    freeType(ddtype_recvfrom_S);
    freeType(ddtype_sendto_S);
    freeType(ddtype_recvfrom_N);
}

void MPISubdomain::makeGhostcellDDType() {
    int sizes[2] = {ny_sub+1, nx_sub+1};
    int subs[2]; int starts[2];

    // --- X direction (i) ---
    subs[0] = ny_sub+1; subs[1] = 1; // ------------------

    // send right face  (i = nx_sub-1)
    starts[0] = 0; starts[1] = nx_sub-1;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_E);
    MPI_Type_commit(&ddtype_sendto_E);
    // recv left ghost (i = 0)
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_W);
    MPI_Type_commit(&ddtype_recvfrom_W);

    // send left face   (i = 1)
    starts[1] = 1;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_W);
    MPI_Type_commit(&ddtype_sendto_W);
    // recv right ghost(i = nx_sub)
    starts[1] = nx_sub;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_E);
    MPI_Type_commit(&ddtype_recvfrom_E);

    // --- Y direction (j) ---
    subs[0] = 1; subs[1] = nx_sub+1; // ------------------

    // send top face    (j = ny_sub-1)
    starts[0]=ny_sub-1; starts[1]=0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_N);
    MPI_Type_commit(&ddtype_sendto_N);
    // recv bottom ghost(j = 0)
    starts[0]=0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_S);
    MPI_Type_commit(&ddtype_recvfrom_S);

    // send bottom face (j = 1)
    starts[0]=1;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_S);
    MPI_Type_commit(&ddtype_sendto_S);
    // recv top ghost   (j = ny_sub)
    starts[0]=ny_sub;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_N);
    MPI_Type_commit(&ddtype_recvfrom_N);
}

void MPISubdomain::ghostcellUpdate(double* theta,
                                   const CartComm1D& cx,
                                   const CartComm1D& cy,
                                   const GlobalParams& /*params*/) {
    MPI_Request reqs[8];
    int r=0;

    // X
    MPI_Isend(theta, 1, ddtype_sendto_E, cx.east_rank, 111, cx.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_W, cx.west_rank, 111, cx.comm, &reqs[r++]);
    MPI_Isend(theta, 1, ddtype_sendto_W, cx.west_rank, 222, cx.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_E, cx.east_rank, 222, cx.comm, &reqs[r++]);

    // Y
    MPI_Isend(theta, 1, ddtype_sendto_N, cy.east_rank, 333, cy.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_S, cy.west_rank, 333, cy.comm, &reqs[r++]);
    MPI_Isend(theta, 1, ddtype_sendto_S, cy.west_rank, 444, cy.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_N, cy.east_rank, 444, cy.comm, &reqs[r++]);

    MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);
}

void MPISubdomain::indices(const GlobalParams& /*params*/, int rankx, int npx, int ranky, int npy) {
    std::fill(theta_x_left_index.begin(), theta_x_left_index.end(), 1);
    std::fill(theta_x_right_index.begin(), theta_x_right_index.end(), 1);
    if (rankx==0)       theta_x_left_index[0]       = 0;
    if (rankx==npx-1)   theta_x_right_index[ny_sub-1] = 0;

    std::fill(theta_y_left_index.begin(), theta_y_left_index.end(), 1);
    std::fill(theta_y_right_index.begin(), theta_y_right_index.end(), 1);
    if (ranky==0)       theta_y_left_index[0]       = 0;
    if (ranky==npy-1)   theta_y_right_index[nx_sub-1] = 0;
}

// x_sub가 값이 저장되는 위치정보인줄 알았는데, 첫번째로 값이 음수가 들어감 
// 이유는 모르겠지만 나중에 다시 함 봐야할 듯 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void MPISubdomain::mesh(const GlobalParams& params,
                        int rankx,int ranky,
                        int npx,int npy) {

    double dx = params.lx / (params.nx - 1);
    for(int i=0; i<=nx_sub;++i) {
        x_sub[i] = (ista - 1 + i - 1) * dx;
        dmx_sub[i] = dx;
    }
    double dy = params.ly / params.ny;
    for(int j=0; j<=ny_sub;++j) {
        y_sub[j] = (jsta - 1 + j) * dy;
        dmy_sub[j] = dy;
    }
    // boundary adjustments omitted
}

void MPISubdomain::initialization(double* theta,
                                  const GlobalParams& params,
                                  int ranky,int npy) {
    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;
    double PI = 3.14159265358979323846;
    for(int j=0;j<=ny_sub;++j)
    for(int i=0;i<=nx_sub;++i) {
        int idx =  j * nx1 + i;
        theta[idx] = (params.theta_cold - params.theta_hot) / params.ly * y_sub[j]
                   + params.theta_hot
                   + sin(4 * PI / params.lx * x_sub[i])
                   * sin(4 * PI / params.ly * y_sub[j]);
    }
}

void MPISubdomain::initialization_debug(double* theta,
                                  const GlobalParams& params,
                                  int myrank) {
    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;
    double PI = 3.14159265358979323846;
    for(int j=0;j<=ny_sub;++j)
    for(int i=0;i<=nx_sub;++i) {
        int idx =  j * nx1 + i;
        theta[idx] = myrank;
    }
}

void MPISubdomain::boundary(const double* theta,
                            const GlobalParams& params,
                            int rankx, int npx,
                            int ranky,int npy) {

    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;

    // ghost cell 에 있는 값을 가져온다.
    int idx;
    // x direction
    for (int j=0; j<ny1; ++j) {
        // i=0
        idx = j*nx1 + 0;
        theta_x_left_sub[j] = theta[idx];

        // i=nx1-1
        idx = j*nx1 + (nx_sub);
        theta_x_right_sub[j] = theta[idx];
    }

    // y
    for (int i=0; i<nx1; ++i) {
        // j=0
        idx = 0*nx1 + i;
        theta_y_left_sub[i] = theta[idx];

        // j = ny1-1
        idx = (ny_sub)*nx1 + i;
        theta_y_right_sub[i] = theta[idx];


    }

    // apply Dirichlet BC at physical walls
    if (ranky==0) {
        for (int i=0; i<nx1; ++i) {
            theta_y_left_sub[i] = params.theta_y_D;
        }
    }
    if (ranky==npy-1) {
        for (int i=0; i<nx1; ++i) {
            theta_y_right_sub[i] = params.theta_y_D;
        }
    }

    if (rankx==0) {
        for (int j=0; j<ny1; ++j) {
                theta_x_left_sub[j] = params.theta_x_D;
        }
    }
    if (rankx==npx-1) {
        for (int j=0; j<ny1; ++j) {
                theta_x_right_sub[j] = params.theta_x_D;
        }
    }
}

MPISubdomain sub;