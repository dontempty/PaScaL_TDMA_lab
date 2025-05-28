// mpi_subdomain.cpp
#include "mpi_subdomain.hpp"
#include <cmath>
#include <algorithm>

static void para_range(int start, int end, int nproc, int rank, int &sta, int &iend) {
    int len = end - start + 1;
    int base = len / nproc;
    int rem  = len % nproc;
    sta = start + rank * base + std::min(rank, rem);
    iend = sta + base - 1 + (rank < rem ? 1 : 0);
}

void MPISubdomain::make(const GlobalParams& params,
                        int npx, int rankx,
                        int npy, int ranky,
                        int npz, int rankz) {
    para_range(0, params.nx-2, npx, rankx, ista, iend);
    nx_sub = iend - ista + 2;
    para_range(0, params.ny-2, npy, ranky, jsta, jend);
    ny_sub = jend - jsta + 2;
    para_range(0, params.nz-2, npz, rankz, ksta, kend);
    nz_sub = kend - ksta + 2;

    x_sub.resize(nx_sub+1);
    dmx_sub.resize(nx_sub+1);
    y_sub.resize(ny_sub+1);
    dmy_sub.resize(ny_sub+1);
    z_sub.resize(nz_sub+1);
    dmz_sub.resize(nz_sub+1);

    thetaBC3_sub.assign((nx_sub+1)*(nz_sub+1),0.0);
    thetaBC4_sub.assign((nx_sub+1)*(nz_sub+1),0.0);
    jmbc_index.assign(ny_sub+1,1);
    jpbc_index.assign(ny_sub+1,1);
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
    freeType(ddtype_sendto_F);
    freeType(ddtype_recvfrom_B);
    freeType(ddtype_sendto_B);
    freeType(ddtype_recvfrom_F);
}

void MPISubdomain::makeGhostcellDDType() {
    int sizes[3] = { nx_sub+1, ny_sub+1, nz_sub+1};
    int subs[3]; int starts[3];

    // --- X direction (i) ---
    subs[0] = 1; subs[1] = ny_sub+1; subs[2] = nz_sub+1;
    // send right face  (i = nx_sub-1)
    starts[0] = nx_sub-1; starts[1] = 0; starts[2] = 0;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_E);
    MPI_Type_commit(&ddtype_sendto_E);
    // recv left ghost (i = 0)
    starts[0] = 0;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_W);
    MPI_Type_commit(&ddtype_recvfrom_W);

    // send left face   (i = 1)
    starts[0] = 1;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_W);
    MPI_Type_commit(&ddtype_sendto_W);
    // recv right ghost(i = nx_sub)
    starts[0] = nx_sub;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_E);
    MPI_Type_commit(&ddtype_recvfrom_E);

    // --- Y direction (j) ---
    subs[0] = nx_sub+1; subs[1] = 1; subs[2] = nz_sub+1;
    // send top face    (j = ny_sub-1)
    starts[0]=0; starts[1]=ny_sub-1; starts[2]=0;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_N);
    MPI_Type_commit(&ddtype_sendto_N);
    // recv bottom ghost(j = 0)
    starts[1]=0;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_S);
    MPI_Type_commit(&ddtype_recvfrom_S);

    // send bottom face (j = 1)
    starts[1]=1;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_S);
    MPI_Type_commit(&ddtype_sendto_S);
    // recv top ghost   (j = ny_sub)
    starts[1]=ny_sub;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_N);
    MPI_Type_commit(&ddtype_recvfrom_N);

    // --- Z direction (k) ---
    subs[0] = nx_sub+1; subs[1] = ny_sub+1; subs[2] = 1;
    // send front face  (k = nz_sub-1)
    starts[0]=0; starts[1]=0; starts[2]=nz_sub-1;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_F);
    MPI_Type_commit(&ddtype_sendto_F);
    // recv back ghost  (k = 0)
    starts[2]=0;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_B);
    MPI_Type_commit(&ddtype_recvfrom_B);

    // send back face   (k = 1)
    starts[2]=1;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_sendto_B);
    MPI_Type_commit(&ddtype_sendto_B);
    // recv front ghost (k = nz_sub)
    starts[2]=nz_sub;
    MPI_Type_create_subarray(3, sizes, subs, starts, MPI_ORDER_C,
                             MPI_DOUBLE, &ddtype_recvfrom_F);
    MPI_Type_commit(&ddtype_recvfrom_F);
}

void MPISubdomain::ghostcellUpdate(double* theta,
                                   const CartComm1D& cx,
                                   const CartComm1D& cy,
                                   const CartComm1D& cz,
                                   const GlobalParams& /*params*/) {
    MPI_Request reqs[12];
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

    // Z
    MPI_Isend(theta, 1, ddtype_sendto_F, cz.east_rank, 555, cz.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_B, cz.west_rank, 555, cz.comm, &reqs[r++]);
    MPI_Isend(theta, 1, ddtype_sendto_B, cz.west_rank, 666, cz.comm, &reqs[r++]);
    MPI_Irecv(theta, 1, ddtype_recvfrom_F, cz.east_rank, 666, cz.comm, &reqs[r++]);

    MPI_Waitall(r, reqs, MPI_STATUSES_IGNORE);
}

void MPISubdomain::indices(const GlobalParams& /*params*/, int ranky, int npy) {
    std::fill(jmbc_index.begin(), jmbc_index.end(), 1);
    std::fill(jpbc_index.begin(), jpbc_index.end(), 1);
    if (ranky==0)       jmbc_index[0]       = 0;
    if (ranky==npy-1)   jpbc_index[ny_sub-1] = 0;
}

void MPISubdomain::mesh(const GlobalParams& params,
                        int rankx,int ranky,int rankz,
                        int npx,int npy,int npz) {
    double dx = params.lx / (params.nx - 1);
    for(int i=0;i<=nx_sub;++i) {
        x_sub[i] = (ista - 1 + i - 1) * dx;
        dmx_sub[i] = dx;
    }
    double dy = params.ly / params.ny;
    for(int j=0;j<=ny_sub;++j) {
        y_sub[j] = (jsta - 1 + j) * dy;
        dmy_sub[j] = dy;
    }
    double dz = params.lz / (params.nz - 1);
    for(int k=0;k<=nz_sub;++k) {
        z_sub[k] = (ksta - 1 + k - 1) * dz;
        dmz_sub[k] = dz;
    }
    // boundary adjustments omitted
}

void MPISubdomain::initialization(double* theta,
                                  const GlobalParams& params,
                                  int ranky,int npy, int myrank) {
    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;
    double PI = 3.14159265358979323846;
    for(int k=0;k<=nz_sub;++k)
    for(int j=0;j<=ny_sub;++j)
    for(int i=0;i<=nx_sub;++i) {
        int idx = (k * ny1 + j) * nx1 + i;
        theta[idx] = myrank;
        // theta[idx] = (params.theta_cold - params.theta_hot) / params.ly * y_sub[j]
        //            + params.theta_hot
        //            + sin(4 * PI / params.lx * x_sub[i])
        //            * sin(4 * PI / params.lz * z_sub[k])
        //            * sin(4 * PI / params.ly * y_sub[j]);
    }
}

void MPISubdomain::boundary(const double* theta,
                            const GlobalParams& params,
                            int ranky,int npy) {
    int nx1 = nx_sub + 1;
    int ny1 = ny_sub + 1;
    // Normal subdomain: extract ghostcell from field
    for(int k=0; k<=nz_sub; ++k) {
        for(int i=0; i<=nx_sub; ++i) {
            // lower wall (j=0)
            int idx_low = (k*ny1 + 0)*nx1 + i;
            thetaBC3_sub[k*nx1 + i] = theta[idx_low];
            // upper wall (j=ny_sub)
            int idx_high = (k*ny1 + ny_sub)*nx1 + i;
            thetaBC4_sub[k*nx1 + i] = theta[idx_high];
        }
    }
    // apply Dirichlet BC at physical walls
    if(ranky == 0) {
        for(int k=0; k<=nz_sub; ++k)
            for(int i=0; i<=nx_sub; ++i)
                thetaBC3_sub[k*nx1 + i] = params.theta_hot;
    }
    if(ranky == npy - 1) {
        for(int k=0; k<=nz_sub; ++k)
            for(int i=0; i<=nx_sub; ++i)
                thetaBC4_sub[k*nx1 + i] = params.theta_cold;
    }
}