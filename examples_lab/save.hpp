// save.hpp
#ifndef SAVE_HPP
#define SAVE_HPP

#include <filesystem>   // C++17
#include <fstream>
#include <vector>
#include <string>
#include <iostream>

void save_rhs_to_csv(const std::vector<double>& rhs,
                     int nx, int ny,
                     const std::string& folder,
                     const std::string& filename)
{
    // 1) 폴더가 없다면 생성
    std::filesystem::path dir(folder);
    if (!std::filesystem::exists(dir)) {
        if (!std::filesystem::create_directories(dir)) {
            std::cerr << "Error: 디렉토리 생성 실패: " << folder << std::endl;
            return;
        }
    }

    // 2) 파일 경로 조합
    std::filesystem::path filepath = dir / filename;

    // 3) 파일 열기
    std::ofstream ofs(filepath.string());
    if (!ofs.is_open()) {
        std::cerr << "Error: 파일을 열 수 없습니다: " << filepath << std::endl;
        return;
    }

    // 4) CSV로 쓰기
    for (int j = 0; j < ny; ++j) {
        for (int i = 0; i < nx; ++i) {
            ofs << rhs[j*nx + i];
            if (i < nx - 1) ofs << ',';
        }
        ofs << '\n';
    }

    ofs.close();
    std::cout << "Saved rhs to " << filepath << std::endl;
}


#endif // SAVE_HPP