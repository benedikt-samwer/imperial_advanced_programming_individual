#include <iostream>
#include <chrono>
#include <cassert>
#include <random>
#include <stdexcept>

// Include your matrix header
#include "Matrix_06031927.h"

// Use the namespace to avoid prefixing every Matrix_06031927 with adv_prog_cw::
using namespace adv_prog_cw;
using namespace std;
using namespace std::chrono;
typedef double fT; // Floating-point type for matrix elements

// Function to generate a random matrix
Matrix_06031927<fT> generateRandomMatrix(size_t rows, size_t cols) {
    Matrix_06031927<fT> mat(rows, cols);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<fT> dis(-10.0, 10.0); // Random values between -10 and 10
    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat(i, j) = dis(gen);
        }
    }
    return mat;
}

// Function to check if two matrices are approximately equal
bool matricesApproximatelyEqual(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b, fT tol = 1e-6) {
    if (a.Rows() != b.Rows() || a.Cols() != b.Cols()) return false;
    for (size_t i = 0; i < a.Rows(); i++) {
        for (size_t j = 0; j < a.Cols(); j++) {
            if (abs(a(i, j) - b(i, j)) > tol) return false;
        }
    }
    return true;
}

// Function to test accuracy of scalar multiplication
void testScalarMultiplicationAccuracy(size_t size) {
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> original = mat;
    fT scalar = 2.5;
    mat *= scalar;
    mat *= (1.0 / scalar); // Reverse the multiplication
    assert(matricesApproximatelyEqual(mat, original));
    cout << "Scalar multiplication accuracy test passed for " << size << "x" << size << " matrix." << endl;
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
        Matrix_06031927<fT> product = mat * inv;
        Matrix_06031927<fT> identity(size, size);
        identity.Identity();
        assert(matricesApproximatelyEqual(product, identity));
        cout << "Matrix inversion accuracy test passed for " << size << "x" << size << " matrix." << endl;
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
    // Test sizes
    vector<size_t> sizes = {10, 100, 1000, 10000};
/*
    // Accuracy Tests
    cout << "\n=== Accuracy Tests ===\n";
    testDeterminantAccuracy(); // Test determinant with a small known matrix
    for (size_t size : sizes) {
        testScalarMultiplicationAccuracy(size);
        testInversionAccuracy(size);
    }
*/
    // Speed Tests with Increasing Size and Complexity
    cout << "\n=== Speed Tests ===\n";
    for (size_t size : sizes) {
        cout << "\nTesting " << size << "x" << size << " matrix:\n";

        // Generate a random matrix for testing
        Matrix_06031927<fT> mat = generateRandomMatrix(size, size);

        // Test Scalar Multiplication Speed
        measureTime([&]() { mat *= 2.5; }, "Scalar multiplication", size);

        // Test Determinant Speed
        measureTime([&]() { mat.Determinant(); }, "Determinant calculation", size);

        // Test Inversion Speed
        Matrix_06031927<fT> inv;
        measureTime([&]() { mat.Inverse(inv); }, "Matrix inversion", size);
    }

    cout << "\nAll tests completed successfully." << endl;
    return 0;
}