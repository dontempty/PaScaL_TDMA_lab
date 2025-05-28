// global.cpp (단순 디버깅 버전)
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <array>
#include <map>
#include "headers.hpp"

// Define global variables
double Pr, Ra;
int Tmax;
int nx, ny, nz;
int nxm, nym, nzm;
int nxp, nyp, nzp;
double dt, dtStart, tStart;
double lx, ly, lz;
double dx, dy, dz;
double theta_cold, theta_hot, alphaG, nu, Ct;

void global_inputpara(const std::string& filename, std::array<int, 3>& np_dim) {
    // std::cout << "=== Starting global_inputpara ===" << std::endl;
    // std::cout << "Filename: " << filename << std::endl;
    
    // 파일 스트림 생성
    std::ifstream infile;
    infile.open(filename.c_str());
    
    if (!infile.is_open()) {
        std::cerr << "ERROR: Cannot open file: " << filename << std::endl;
        throw std::runtime_error("Input file could not be opened: " + filename);
    }
    
    // std::cout << "File opened successfully!" << std::endl;

    // std::map 사용 (더 안전)
    std::map<std::string, std::string> param_map;
    std::string line;
    int line_count = 0;
    
    while (std::getline(infile, line)) {
        line_count++;
        
        if (line.empty()) continue;

        // 주석 제거
        size_t exclam = line.find('!');
        if (exclam != std::string::npos) {
            line = line.substr(0, exclam);
        }

        std::istringstream iss(line);
        std::string key, value;
        if (iss >> key >> value) {
            if (!key.empty() && !value.empty()) {
                param_map[key] = value;
                // std::cout << "Parsed: " << key << " = " << value << std::endl;
            }
        }
    }
    
    infile.close();
    // std::cout << "File closed. Total lines processed: " << line_count << std::endl;

    // 필수 파라미터 확인 및 할당
    try {
        nx = std::stoi(param_map.at("nx"));
        ny = std::stoi(param_map.at("ny"));
        nz = std::stoi(param_map.at("nz"));
        
        np_dim[0] = std::stoi(param_map.at("npx"));
        np_dim[1] = std::stoi(param_map.at("npy"));
        np_dim[2] = std::stoi(param_map.at("npz"));
        
        Tmax = std::stoi(param_map.at("Tmax"));
        
    } catch (const std::out_of_range& e) {
        std::cerr << "Missing required parameter: " << e.what() << std::endl;
        throw std::runtime_error("Missing required parameter");
    } catch (const std::invalid_argument& e) {
        std::cerr << "Invalid parameter value: " << e.what() << std::endl;
        throw std::runtime_error("Invalid parameter value");
    }

    // 나머지 초기화
    Pr = 5.0; 
    Ra = 200.0;
    nx += 1; ny += 1; nz += 1;
    nxm = nx - 1; nym = ny - 1; nzm = nz - 1;
    nxp = nx + 1; nyp = ny + 1; nzp = nz + 1;
    dtStart = 5.0e-3; 
    tStart = 0.0;
    lx = 1.0; ly = 1.0; lz = 1.0;
    theta_cold = -1.0;
    theta_hot = 2.0 + theta_cold;
    alphaG = 1.0;
    
    double denominator = alphaG * Pr * std::pow(ly, 3.0) * (theta_hot - theta_cold);
    // std::cout << "Denominator for nu calculation: " << denominator << std::endl;
    
    nu = 1.0 / std::sqrt(Ra / denominator);
    Ct = nu / Pr;
    
    // std::cout << "=== Parameter parsing completed successfully! ===" << std::endl;
}