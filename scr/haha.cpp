// #include <iostream>
// #include <vector>
// #include <chrono>
// #include <numeric>  // for std::iota if needed

// int main() {
//     const size_t nprocs = 1000;

//     // 방법 1: resize + fill
//     std::vector<int> vec1;
//     auto t1_start = std::chrono::high_resolution_clock::now();
//     vec1.resize(nprocs);
//     std::fill(vec1.begin(), vec1.end(), 1);
//     auto t1_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed1 = t1_end - t1_start;

//     // 방법 2: 복사 대입 (deep copy처럼)
//     auto t2_start = std::chrono::high_resolution_clock::now();
//     std::vector<int> vec2(nprocs, 1);
//     auto t2_end = std::chrono::high_resolution_clock::now();
//     std::chrono::duration<double> elapsed2 = t2_end - t2_start;

//     // 출력
//     std::cout << "Method 1 (resize + fill): " << elapsed1.count() << " seconds\n";
//     std::cout << "Method 2 (copy-init):     " << elapsed2.count() << " seconds\n";

//     return 0;
// }

#include <iostream>
#include <vector>
#include <chrono>

using namespace std;
using namespace std::chrono;

void test_resize(int outer, int inner, int n_repeat) {
    vector<vector<double>> mat;
    auto start = high_resolution_clock::now();

    for (int r = 0; r < n_repeat; ++r) {
        mat.resize(outer);
        for (int i = 0; i < outer; ++i)
            mat[i].resize(inner);
    }

    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    cout << "Resize time: " << elapsed.count() << " seconds\n";
}

void test_redeclare(int outer, int inner, int n_repeat) {
    vector<vector<double>> mat;
    auto start = high_resolution_clock::now();

    for (int r = 0; r < n_repeat; ++r) {
        mat = vector<vector<double>>(outer, vector<double>(inner));
    }

    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    cout << "Re-declare time: " << elapsed.count() << " seconds\n";
}

int main() {
    int outer = 10;
    int inner = 10;
    int n_repeat = 10;

    test_resize(outer, inner, n_repeat);
    test_redeclare(outer, inner, n_repeat);

    return 0;
}

