// global.cpp 에서 load() 부분만 수정

#include "global.hpp"
#include <fstream>
#include <sstream>
#include <map>
#include <cmath>
#include <algorithm>
#include <cctype>

// trim 유틸
static inline void trim(std::string &s) {
    auto f = [](char c){ return !std::isspace((unsigned char)c); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), f));
    s.erase(std::find_if(s.rbegin(), s.rend(), f).base(), s.end());
}

void GlobalParams::load(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) throw std::runtime_error("Cannot open input file: " + filename);

    std::map<std::string,std::string> param;
    std::string line;
    while (std::getline(infile, line)) {
        // 1) 주석 제거
        if (auto pos = line.find('!'); pos != std::string::npos)
            line.erase(pos);

        // 2) 앞뒤 공백 제거
        trim(line);
        if (line.empty()) continue;

        std::string key, value;
        auto eq = line.find('=');
        if (eq != std::string::npos) {
            // key = value 형식
            key   = line.substr(0, eq);
            value = line.substr(eq + 1);
        } else {
            // key value 형식
            std::istringstream iss(line);
            if (!(iss >> key >> value)) continue;
        }
        trim(key); trim(value);
        if (key.empty() || value.empty()) continue;
        param[key] = value;
    }

    // 이제 모든 형식이 param에 들어가 있음
    nx       = std::stoi(param.at("nx"));
    ny       = std::stoi(param.at("ny"));
    // nz       = std::stoi(param.at("nz"));
    np_dim[0]= std::stoi(param.at("npx"));
    np_dim[1]= std::stoi(param.at("npy"));
    // np_dim[2]= std::stoi(param.at("npz"));
    Tmax     = std::stoi(param.at("Tmax"));

#ifdef THREAD_PARAMS
    thread_in_x        = std::stoi(param.at("thread_in_x"));
    thread_in_y        = std::stoi(param.at("thread_in_y"));
    thread_in_z        = std::stoi(param.at("thread_in_z"));
    thread_in_x_pascal = std::stoi(param.at("thread_in_x_pascal"));
    thread_in_y_pascal = std::stoi(param.at("thread_in_y_pascal"));
#endif

    // 고정값 초기화
    Pr = 5.0; Ra = 200.0;
    nx++; ny++; // nz++;   // 1 더해준다.
    nxm = nx - 1; nym = ny - 1; // nzm = nz - 1;
    nxp = nx + 1; nyp = ny + 1; // nzp = nz + 1;

    dtStart = 0.01; tStart = 0.0;
    x0 = -1; xN = 1;
    y0 = -1; yN = 1;
    z0 = -1; zN = 1;
    lx = (xN - x0); ly = (yN - y0); lz = (zN - z0);

    // theta_cold = -1.0;
    // theta_hot  =  2.0 + theta_cold;
    theta_x_L_D = 0;
    theta_x_R_D = 0;
    theta_y_L_D = 0;
    theta_y_R_D = 0;

    alphaG     =  1.0;

    double denom = alphaG * Pr * std::pow(ly, 3) * (theta_hot - theta_cold);
    nu = 1.0 / std::sqrt(Ra / denom);
    Ct = nu / Pr;

    // dx, dy, dz 계산 (옵션)
    dx = lx / (nx - 1);
    dy = ly / (ny - 1);
    // dz = lz / (nz - 1);
}

GlobalParams params;
