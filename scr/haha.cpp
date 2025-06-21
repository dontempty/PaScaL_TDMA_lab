#include <iostream>
#include <vector>
#include <chrono>
#include <numeric>  // for std::iota if needed

int main() {
    const size_t nprocs = 1000;

    // 방법 1: resize + fill
    std::vector<int> vec1;
    auto t1_start = std::chrono::high_resolution_clock::now();
    vec1.resize(nprocs);
    std::fill(vec1.begin(), vec1.end(), 1);
    auto t1_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed1 = t1_end - t1_start;

    // 방법 2: 복사 대입 (deep copy처럼)
    auto t2_start = std::chrono::high_resolution_clock::now();
    std::vector<int> vec2(nprocs, 1);
    auto t2_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = t2_end - t2_start;

    // 출력
    std::cout << "Method 1 (resize + fill): " << elapsed1.count() << " seconds\n";
    std::cout << "Method 2 (copy-init):     " << elapsed2.count() << " seconds\n";

    return 0;
}
