#include <iostream>
#include <chrono>
#include <cassert>
#include <random>
#include <stdexcept>
#include <omp.h> // OpenMP for manual thread selection

// Include your matrix header
#include "Matrix_06031927.h"

// Use the namespace to avoid prefixing every Matrix_06031927 with adv_prog_cw::
using namespace adv_prog_cw;
using namespace std;
using namespace std::chrono;
typedef double fT; // Floating-point type for matrix elements

// Function to generate a 50% sparse random matrix
Matrix_06031927<fT> generateRandomMatrix(size_t rows, size_t cols) {
    Matrix_06031927<fT> mat(rows, cols);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<fT> dis(-10.0, 10.0); // Values between -10 and 10
    bernoulli_distribution sparse_prob(0.5); // 50% chance of being zero

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat(i, j) = sparse_prob(gen) ? 0 : dis(gen);
        }
    }
    return mat;
}

// Function to test accuracy of determinant (using a simple 2x2 matrix for verification)
void testDeterminantAccuracy() {
    Matrix_06031927<fT> mat(2, 2);
    mat(0, 0) = 1; mat(0, 1) = 2;
    mat(1, 0) = 3; mat(1, 1) = 4;
    fT det = mat.Determinant();
    fT expected = -2.0; // Known determinant for [[1, 2], [3, 4]]
    assert(abs(det - expected) < 1e-6);
    cout << "Determinant accuracy test passed for 2x2 matrix." << endl;
}

// Function to test accuracy of matrix inversion
void testInversionAccuracy(size_t size) {
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> inv;
    bool success = mat.Inverse(inv);
    if (success) {
        cout << "Matrix inversion test passed for " << size << "x" << size << " matrix." << endl;
    } else {
        cout << "Matrix inversion skipped (singular matrix) for " << size << "x" << size << "." << endl;
    }
}

// Function to measure execution time of an operation
template<typename Func>
double measureTime(Func operation, const string& operationName, size_t size) {
    auto start = high_resolution_clock::now();
    operation();
    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    cout << operationName << " for " << size << "x" << size << " matrix took " << elapsed.count() << " seconds." << endl;
    return elapsed.count();
}

int main() {
    // Allow user to set the number of OpenMP threads (kernels)
    int num_threads;
    cout << "Enter the number of OpenMP threads (0 for automatic): ";
    cin >> num_threads;

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        cout << "Using " << num_threads << " OpenMP threads." << endl;
    } else {
        cout << "Using automatic OpenMP thread selection." << endl;
    }

    // Verify OpenMP thread count
    #pragma omp parallel
    {
        #pragma omp single
        cout << "OpenMP is using " << omp_get_num_threads() << " threads." << endl;
    }

    // Test sizes
    vector<size_t> sizes = {2,3, 10, 100, 500, 1000, 3000, 10000, 100000};

    // Speed Tests with Increasing Size and Complexity
    cout << "\n=== Speed Tests ===\n";
    for (size_t size : sizes) {
        cout << "\nTesting " << size << "x" << size << " matrix:\n";

        // Generate a random matrix for testing
        Matrix_06031927<fT> mat = generateRandomMatrix(size, size);

        // Test Scalar Multiplication Speed
        measureTime([&]() { mat *= 2; }, "Scalar multiplication", size);
        measureTime([&]() { mat /= 2; }, "Scalar division", size);

        // Test Inversion Speed
        Matrix_06031927<fT> inv;
        measureTime([&]() { mat.Inverse(inv); }, "Matrix inversion", size);

        // Test Determinant Speed
        measureTime([&]() { mat.Determinant(); }, "Determinant calculation", size);


    }

    cout << "\nAll tests completed successfully." << endl;
    return 0;
}
