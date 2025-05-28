// headers.hpp
#ifndef HEADERS_HPP
#define HEADERS_HPP

#include <string>
#include <vector>
#include <array>

// headers.hpp

// Global values from global.cpp
extern double Pr, Ra, nu, Ct;
extern int Tmax;
extern int nx, ny, nz;
extern int nxm, nym, nzm;
extern int nxp, nyp, nzp;
extern double dt, dtStart, tStart;
extern double lx, ly, lz;
extern double dx, dy, dz;
extern double theta_cold, theta_hot, alphaG;

// Function declarations
void global_inputpara(const std::string& filename, std::array<int, 3>& np_dim);


#endif // HEADERS_HPP