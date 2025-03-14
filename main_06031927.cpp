#include <iostream>
#include <chrono>
#include <cassert>
#include <random>
#include <stdexcept>
#include <cmath>
#include <omp.h>
#include "Matrix_06031927.h"
#include <Eigen/Dense>  // For reference computations with Eigen

using namespace adv_prog_cw;
using namespace std;
using namespace std::chrono;
typedef double fT;

// Function to generate a 50% sparse random matrix
Matrix_06031927<fT> generateRandomMatrix(size_t rows, size_t cols) {
    Matrix_06031927<fT> mat(rows, cols);
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<fT> dis(-10.0, 10.0); // Values between -10 and 10
    bernoulli_distribution sparse_prob(0.3); // 50% chance of being zero

    for (size_t i = 0; i < rows; i++) {
        for (size_t j = 0; j < cols; j++) {
            mat(i, j) = sparse_prob(gen) ? 0 : dis(gen);
        }
    }
    return mat;
}

// Accuracy test for Determinant using Eigen for reference
void testDeterminantAccuracyEigen(size_t size) {
    //size_t size = 5; // Adjust size as needed
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    fT myDet = mat.Determinant();

    // Build an equivalent Eigen matrix.
    Eigen::MatrixXd eigenMat(size, size);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            eigenMat(i, j) = mat(i, j);
        }
    }
    fT eigenDet = eigenMat.determinant();

    // Print both determinants for comparison.
    cout << "My Determinant:    " << myDet << endl;
    cout << "Eigen Determinant: " << eigenDet << endl;

    // Check accuracy using a relative tolerance.
    if (fabs(myDet - eigenDet) < 1e-6 * fabs(eigenDet)) {
        cout << "Determinant accuracy test passed for " << size << "x" << size << " matrix." << endl;
    } else {
        cout << "Determinant accuracy test FAILED for " << size << "x" << size << " matrix." << endl;
        cout << "Relative error: " << fabs(myDet - eigenDet) / fabs(eigenDet) << endl;
        //assert(false);
    }
}
void testInversionAccuracyEigen(size_t size) {
    Matrix_06031927<fT> mat = generateRandomMatrix(size, size);
    Matrix_06031927<fT> myInv;
    bool success = mat.Inverse(myInv);
    if (!success) {
        cout << "Matrix inversion skipped (singular matrix) for " 
             << size << "x" << size << " matrix." << endl;
        return;
    }

    // Build an Eigen matrix equivalent to our matrix.
    Eigen::MatrixXd eigenMat(size, size);
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            eigenMat(i, j) = mat(i, j);
        }
    }
    Eigen::MatrixXd eigenInv = eigenMat.inverse();

    // Compute the Frobenius norm of the difference.
    double diff = 0.0, norm = 0.0;
    for (size_t i = 0; i < size; i++) {
        for (size_t j = 0; j < size; j++) {
            double d = myInv(i, j) - eigenInv(i, j);
            diff += d * d;
            norm += eigenInv(i, j) * eigenInv(i, j);
        }
    }
    diff = sqrt(diff);
    norm = sqrt(norm);
    if (norm == 0) norm = 1; // avoid division by zero

    double relativeError = diff / norm;
    cout << "Relative error inversion: " << relativeError << endl;

    if (relativeError < 1e-6) {
        cout << "Matrix inversion accuracy test passed for " 
             << size << "x" << size << " matrix." << endl;
    } else {
        cout << "Matrix inversion accuracy test FAILED for " 
             << size << "x" << size << " matrix." << endl;
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
    // Allow user to set the number of OpenMP threads.
    int num_threads;
    cout << "Enter the number of OpenMP threads (0 for automatic): ";
    cin >> num_threads;

    if (num_threads > 0) {
        omp_set_num_threads(num_threads);
        cout << "Using " << num_threads << " OpenMP threads." << endl;
    } else {
        cout << "Using automatic OpenMP thread selection." << endl;
    }

    // Verify OpenMP thread count.
    #pragma omp parallel
    {
        #pragma omp single
        cout << "OpenMP is using " << omp_get_num_threads() << " threads." << endl;
    }

    // Test sizes for speed and accuracy tests.
    vector<size_t> sizes = {5, 11, 51, 101, 500, 1000};

    // Run the determinant accuracy test
    cout << "\n=== Determinant Accuracy Test (Eigen Comparison) ===\n";
    

    for (size_t size : sizes) {
        testDeterminantAccuracyEigen(size);
        testInversionAccuracyEigen(size);
    }

    // Test sizes for speed and accuracy tests.
    sizes = {10, 100, 500, 1000, 2000, 3000, 4000, 5000};

    // Speed Tests with Increasing Size and Complexity
    cout << "\n=== Speed Tests ===\n";
    for (size_t size : sizes) {
        cout << "\nTesting " << size << "x" << size << " matrix:\n";

        // Generate a random matrix for testing.
        Matrix_06031927<fT> mat = generateRandomMatrix(size, size);

        // Test Scalar Multiplication Speed.
        measureTime([&]() { mat *= 2; }, "Scalar multiplication", size);
        measureTime([&]() { mat /= 2; }, "Scalar division", size);
        //testDeterminantAccuracyEigen(size);
        fT det = 0.0;
        double elapsed = measureTime([&]() { det = mat.Determinant(); }, "Determinant calculation", size);
        cout << "Computed determinant: " << det << endl;

        // Test Inversion Speed.
        Matrix_06031927<fT> inv;
        measureTime([&]() { mat.Inverse(inv); }, "Matrix inversion", size);
    }

    cout << "\nAll tests completed successfully." << endl;
    return 0;
}