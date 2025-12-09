# Advanced Programming Assessment

**Imperial College London**

**MSc Applied Computational Science and Engineering (ACSE)**

**Module:** Advanced Programming (2024/25)

## Overview

This repository contains the individual coursework for the Advanced Programming module. The project focuses on implementing a high-performance, template-based dense matrix class in C++ (`Matrix_06031927`) optimized for speed and numerical stability.

The core objective was to extend a basic matrix class skeleton with advanced linear algebra operations—specifically scalar multiplication/division, determinant calculation, and matrix inversion—while ensuring the implementation is efficient enough to handle large matrices (up to 5000x5000) effectively.

## Key Features & Implementation Details

### 1. Template-Based Design
The `Matrix_06031927` class is templated to support various floating-point types (e.g., `double`, `float`), ensuring flexibility and type safety. It utilizes `std::vector<std::vector<T>>` for dynamic memory management.

### 2. High-Performance Algorithms
The solution leverages **LU Decomposition with Partial Pivoting** as the backbone for both determinant calculation and matrix inversion. This approach was chosen for its $O(n^3)$ complexity and numerical stability compared to other methods like Gaussian elimination without pivoting.

#### Optimizations Implemented:
*   **OpenMP Parallelization**: 
    *   **Scalar Operations**: Element-wise multiplication and division are parallelized using `#pragma omp parallel for collapse(2)` with a threshold check (operations on matrices smaller than $100 \times 100$ run serially to avoid thread overhead).
    *   **LU Decomposition**: The pivot search and the update of the trailing submatrix are parallelized.
    *   **Matrix Inversion**: The solver computes inverse columns in parallel (`#pragma omp parallel for`), treating each column solution ($Ax = e_i$) as an independent task.
    
*   **Vectorization (SIMD)**: 
    *   Inner loops in the LU decomposition and substitution steps are annotated with `#pragma omp simd` to allow the compiler to generate vector instructions (AVX/SSE), significantly speeding up floating-point arithmetic.

*   **Cache Optimization (Blocking)**:
    *   The `ParallelLU` method implements **tiled (blocked) processing** (block size = 128). This technique improves cache locality by keeping working datasets within the CPU cache, reducing expensive memory access latency during the decomposition of large matrices.

*   **Numerical Stability**:
    *   **Partial Pivoting**: Swaps rows to place the largest element in the current column on the diagonal, minimizing numerical error.
    *   **Log-Determinant**: The determinant is computed using the sum of logarithms of the diagonal elements of $U$. This prevents arithmetic overflow/underflow when dealing with very large or very small determinant values.

### 3. Verification & Testing
The project includes a robust testing framework (`main_06031927.cpp`) that:
*   **Validates Correctness**: Compares the results of the custom implementation against the industry-standard **Eigen** library to ensure high precision (checking relative error < $10^{-6}$).
*   **Benchmarks Performance**: Measures execution time for operations on matrices of increasing size (from $5 \times 5$ to $5000 \times 5000$).
*   **Scalability**: Demonstrates the effectiveness of parallelization by allowing the user to configure the number of OpenMP threads.

## Build and Run

To compile the project, ensure you have a C++ compiler supporting C++11 (or later) and OpenMP. You will also need the Eigen library headers in your include path for the tests.

### Compilation
```bash
g++ -O3 -fopenmp -I/path/to/eigen main_06031927.cpp Matrix_06031927.cpp -o matrix_solver
```

### Execution
```bash
./matrix_solver
```
The program will prompt you to enter the number of threads (enter `0` for automatic selection) and then proceed to run accuracy and speed tests.

## File Structure

*   `Matrix_06031927.h`: Header file declaring the matrix class template and interface.
*   `Matrix_06031927.cpp`: Implementation of the matrix operations (scalar ops, LU decomposition, determinant, inverse).
*   `main_06031927.cpp`: Main entry point containing test suites and benchmarking logic.

## License

This project is part of an academic assessment for Imperial College London.

