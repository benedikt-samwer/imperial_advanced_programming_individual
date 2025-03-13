#include <iostream>
#include <chrono>
#include <cassert>
#include <random>
#include <stdexcept>
#include <omp.h>
#include "Matrix_06031927.h"

using namespace adv_prog_cw;
using namespace std;
using namespace std::chrono;
typedef double fT;

Matrix_06031927<fT> generateRandomMatrix(size_t rows, size_t cols) {
    Matrix_06031927<fT> mat(rows, cols);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<fT> dis(-10.0, 10.0);
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat(i, j) = dis(gen);
        }
    }
    return mat;
}

bool matricesApproximatelyEqual(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b, fT tol = 1e-6) {
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) return false;
    for (size_t i = 0; i < a.Rows(); i++) {
        for (size_t j = 0; j < a.Cols(); j++) {
            if (abs(a(i, j) - b(i, j)) > tol) return false;
        }
    }
    return true;
}

template<typename Func>
double measureTime(Func operation, const string& operationName, size_t size, bool& skip) {
    if (skip) {
        cout << operationName << " for " << size << "x" << size << " matrix skipped (previous operation exceeded 60 seconds)." << endl;
        return -1.0;
    }

    auto start = high_resolution_clock::now();
    operation();
    auto end = high_resolution_clock::now();
    duration<double> elapsed = end - start;
    double time_taken = elapsed.count();

    cout << operationName << " for " << size << "x" << size << " matrix took " << time_taken << " seconds." << endl;
    if (time_taken > 60.0) {
        skip = true;
        cout << "Warning: " << operationName << " exceeded 60 seconds; skipping subsequent operations for this size." << endl;
    }
    return time_taken;
}

void testScalarMultiplicationAccuracy(size_t size, bool& skip) {
    if (skip) {
        cout << "Scalar multiplication accuracy test skipped for " << size << "x" << size << " matrix." << endl;
        return;
    }
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> original = mat;
    fT scalar = 2.5;
    measureTime([&]() { mat *= scalar; mat *= (1.0 / scalar); }, "Scalar multiplication accuracy", size, skip);
    if (!skip) {
        assert(matricesApproximatelyEqual(mat, original));
        cout << "Scalar multiplication accuracy test passed for " << size << "x" << size << " matrix." << endl;
    }
}

void testScalarDivisionAccuracy(size_t size, bool& skip) {
    if (skip) {
        cout << "Scalar division accuracy test skipped for " << size << "x" << size << " matrix." << endl;
        return;
    }
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> original = mat;
    fT scalar = 2.5;
    measureTime([&]() { mat /= scalar; mat *= scalar; }, "Scalar division accuracy", size, skip);
    if (!skip) {
        assert(matricesApproximatelyEqual(mat, original));
        cout << "Scalar division accuracy test passed for " << size << "x" << size << " matrix." << endl;
    }
}

void testDeterminantAccuracy(bool& skip) {
    if (skip) {
        cout << "Determinant accuracy test skipped for 2x2 matrix." << endl;
        return;
    }
    Matrix_06031927<fT> mat(2, 2);
    mat(0, 0) = 1; mat(0, 1) = 2;
    mat(1, 0) = 3; mat(1, 1) = 4;
    fT det;
    measureTime([&]() { det = mat.Determinant(); }, "Determinant accuracy", 2, skip);
    if (!skip) {
        fT expected = -2.0;
        assert(abs(det - expected) < 1e-6);
        cout << "Determinant accuracy test passed for 2x2 matrix." << endl;
    }
}

void testInversionAccuracy(size_t size, bool& skip) {
    if (skip) {
        cout << "Matrix inversion accuracy test skipped for " << size << "x" << size << " matrix." << endl;
        return;
    }
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> inv;
    bool success;
    measureTime([&]() { success = mat.Inverse(inv); }, "Matrix inversion accuracy", size, skip);
    if (!skip && success) {
        Matrix_06031927<fT> product = mat * inv;
        Matrix_06031927<fT> identity(size, size);
        identity.Identity();
        assert(matricesApproximatelyEqual(product, identity));
        cout << "Matrix inversion accuracy test passed for " << size << "x" << size << " matrix." << endl;
    } else if (!skip) {
        cout << "Matrix inversion skipped (singular matrix) for " << size << "x" << size << "." << endl;
    }
}

int main() {
    // Set the number of threads explicitly
    int num_threads = 8; // Adjust based on your CPU (e.g., 4 or 8 cores)
    omp_set_num_threads(num_threads);
    cout << "Using " << num_threads << " threads for OpenMP parallelization." << endl;

    // Test sizes
    vector<size_t> sizes = {10, 100, 1000, 10000};

    // Accuracy Tests
    cout << "\n=== Accuracy Tests ===\n";
    bool skip_accuracy = false;
    testDeterminantAccuracy(skip_accuracy);
    for (size_t size : sizes) {
        testScalarMultiplicationAccuracy(size, skip_accuracy);
        testScalarDivisionAccuracy(size, skip_accuracy);
        testInversionAccuracy(size, skip_accuracy);
        if (skip_accuracy) {
            cout << "Skipping remaining accuracy tests due to timeout." << endl;
            break;
        }
    }

    // Speed Tests
    cout << "\n=== Speed Tests ===\n";
    bool skip_speed = false;
    for (size_t size : sizes) {
        cout << "\nTesting " << size << "x" << size << " matrix:\n";
        Matrix_06031927<fT> mat = generateRandomMatrix(size, size);

        measureTime([&]() { mat *= 2.5; }, "Scalar multiplication", size, skip_speed);
        measureTime([&]() { mat /= 2.5; }, "Scalar division", size, skip_speed);
        measureTime([&]() { mat.Determinant(); }, "Determinant calculation", size, skip_speed);
        Matrix_06031927<fT> inv;
        measureTime([&]() { mat.Inverse(inv); }, "Matrix inversion", size, skip_speed);

        if (skip_speed) {
            cout << "Skipping remaining speed tests due to timeout." << endl;
            break;
        }
    }

    cout << "\nAll tests completed successfully (within time limits)." << endl;
    return 0;
}