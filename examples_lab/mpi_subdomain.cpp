// mpi_subdomain.cpp
#include "mpi_subdomain.hpp"
#include "../scr/para_range.hpp"
#include <cmath>
#include "iostream"
#include <algorithm>

const double Pi = 3.14159265358979323846;

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

    // 이놈은 x_left는 x 축의 왼쪽, 즉 y축 크기의 bdy정보를 저장
    theta_x_left_sub.assign((ny_sub+1), 0.0);
    theta_x_right_sub.assign((ny_sub+1), 0.0);
    theta_y_left_sub.assign((nx_sub+1), 0.0);
    theta_y_right_sub.assign((nx_sub+1), 0.0);

    // 이놈은 각 축에서 어느 위치가 bdy인지 나타내는 index
    theta_x_left_index.assign(nx_sub+1, 0);
    theta_x_right_index.assign(nx_sub+1, 0);
    theta_y_left_index.assign(ny_sub+1, 0);
    theta_y_right_index.assign(ny_sub+1, 0);
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
    subs[0] = ny_sub+1;
    subs[1] = 1;

    // send right face  (i = nx_sub-1)
    starts[0] = 0;
    starts[1] = nx_sub-1;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_E);
    MPI_Type_commit(&ddtype_sendto_E);
    // recv left ghost (i = 0)
    starts[0] = 0;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_W);
    MPI_Type_commit(&ddtype_recvfrom_W);

    // send left face   (i = 1)
    starts[0] = 0;
    starts[1] = 1;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_W);
    MPI_Type_commit(&ddtype_sendto_W);
    // recv right ghost(i = nx_sub)
    starts[0] = 0;
    starts[1] = nx_sub;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_E);
    MPI_Type_commit(&ddtype_recvfrom_E);

    // --- Y direction (j) ---
    subs[0] = 1;
    subs[1] = nx_sub+1;

    // send top face    (j = ny_sub-1)
    starts[0] = ny_sub-1;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_N);
    MPI_Type_commit(&ddtype_sendto_N);
    // recv bottom ghost(j = 0)
    starts[0] = 0;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_S);
    MPI_Type_commit(&ddtype_recvfrom_S);

    // send bottom face (j = 1)
    starts[0] = 1;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_S);
    MPI_Type_commit(&ddtype_sendto_S);
    // recv top ghost   (j = ny_sub)
    starts[0] = ny_sub;
    starts[1] = 0;
    MPI_Type_create_subarray(2, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_N);
    MPI_Type_commit(&ddtype_recvfrom_N);
}

void MPISubdomain::ghostcellUpdate(std::vector<double>& theta,
                                   const CartComm1D& cx,
                                   const CartComm1D& cy,
                                   const GlobalParams& /*params*/) {
    MPI_Request reqs[8];
    int r=0;

    // X
    MPI_Isend(theta.data(), 1, ddtype_sendto_E, cx.east_rank, 111, cx.comm, &reqs[r++]);
    MPI_Irecv(theta.data(), 1, ddtype_recvfrom_W, cx.west_rank, 111, cx.comm, &reqs[r++]);
    MPI_Isend(theta.data(), 1, ddtype_sendto_W, cx.west_rank, 222, cx.comm, &reqs[r++]);
    MPI_Irecv(theta.data(), 1, ddtype_recvfrom_E, cx.east_rank, 222, cx.comm, &reqs[r++]);

    // Y
    MPI_Isend(theta.data(), 1, ddtype_sendto_N, cy.east_rank, 333, cy.comm, &reqs[r++]);
    MPI_Irecv(theta.data(), 1, ddtype_recvfrom_S, cy.west_rank, 333, cy.comm, &reqs[r++]);
    MPI_Isend(theta.data(), 1, ddtype_sendto_S, cy.west_rank, 444, cy.comm, &reqs[r++]);
    MPI_Irecv(theta.data(), 1, ddtype_recvfrom_N, cy.east_rank, 444, cy.comm, &reqs[r++]);

    MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);
}

void MPISubdomain::indices(const GlobalParams& /*params*/, int rankx, int npx, int ranky, int npy) {
    std::fill(theta_x_left_index.begin(), theta_x_left_index.end(), 0);
    std::fill(theta_x_right_index.begin(), theta_x_right_index.end(), 0);
    if (rankx==0)       theta_x_left_index[1]       = 1;
    if (rankx==npx-1)   theta_x_right_index[nx_sub-1] = 1;

    std::fill(theta_y_left_index.begin(), theta_y_left_index.end(), 0);
    std::fill(theta_y_right_index.begin(), theta_y_right_index.end(), 0);
    if (ranky==0)       theta_y_left_index[1]       = 1;
    if (ranky==npy-1)   theta_y_right_index[ny_sub-1] = 1;
}

void MPISubdomain::mesh(const GlobalParams& params,
                        int rankx,int ranky,
                        int npx,int npy) {

    double dx = params.lx / (params.nx - 1);
    for(int i=0; i<=nx_sub;++i) {
        
        if (rankx==0 && i==0) {
             x_sub[i] = params.x0;
             dmx_sub[i] = dx/2;
        }
        else if (rankx==npx-1 && i==nx_sub) {
            x_sub[i] = params.xN;
            dmx_sub[i] = dx/2;
        }
        else {
            x_sub[i] = params.x0 + dx/2 + (ista - 2 + i)*dx;
            dmx_sub[i] = dx;
        }
    }

    double dy = params.ly / (params.ny-1);
    for(int j=0; j<=ny_sub;++j) {

        if (ranky==0 && j==0) {
             y_sub[j] = params.y0;
             dmy_sub[j] = dy/2;
        }
        else if (ranky==npy-1 && j==ny_sub) {
            y_sub[j] = params.yN;
            dmy_sub[j] = dy/2;
        }
        else {
            y_sub[j] = params.y0 + dy/2 + (jsta - 2 + j)*dy;
            dmy_sub[j] = dy;
        }
    }
    // boundary adjustments omitted
}

void MPISubdomain::initialization(std::vector<double>& theta,
                                  const GlobalParams& params) {
    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;
    int myrank;
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    
    for(int j=0; j<ny1; ++j) {
        for(int i=0; i<nx1; ++i) {
            int idx =  j * nx1 + i;

            theta[idx] = sin(Pi*x_sub[i]) * sin(Pi*y_sub[j]) * exp(-2.0 * Pi*Pi * 0.0) + cos(Pi*x_sub[i]) * cos(Pi*y_sub[j]);
        }
    }
}

void MPISubdomain::boundary(std::vector<double>& theta,
                            const GlobalParams& params,
                            int rankx, int npx,
                            int ranky,int npy) {

    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;

    // ghost cell 에 있는 값을 가져온다.
    int idx;
    // y 축 bdy 정보 담기
    for (int j=0; j<ny1; ++j) {
        // i=0
        idx = j*nx1 + 0;
        theta_x_left_sub[j] = theta[idx];

        // i = nx1-1
        idx = j*nx1 + (nx1-1);
        theta_x_right_sub[j] = theta[idx];
    }
    // x 축 bdy 정보 담기
    for (int i=0; i<nx1; ++i) {
        // j=0
        idx = 0*nx1 + i;
        theta_y_left_sub[i] = theta[idx];

        // j=ny-1
        idx = (ny1-1)*nx1 + i;
        theta_y_right_sub[i] = theta[idx];
    }

    // apply Dirichlet BC at physical walls

    // y 축 bdy 정보 담기
    if (rankx==0) {
        for (int j=0; j<ny1; ++j) {
            theta_x_left_sub[j] = params.theta_y_L_D;
        }
    }
    if (rankx==npx-1) {
        for (int j=0; j<ny1; ++j) {
            theta_x_right_sub[j] = params.theta_y_R_D;
        }
    }

    // x 축 bdy 정보 담기
    if (ranky==0) {
        for (int i=0; i<nx1; ++i) {
                theta_y_left_sub[i] = params.theta_x_L_D;
        }
    }
    if (ranky==npy-1) {
        for (int i=0; i<nx1; ++i) {
            theta_y_right_sub[i] = params.theta_x_R_D;
        }
    }
}

// MPISubdomain sub;