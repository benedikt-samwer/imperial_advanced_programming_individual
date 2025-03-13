// Step 1: RENAME THE _ADV_PROG_Matrix_06031927_H_ Label To _ADV_PROG_Matrix_06031927_00112233_H_
// (Replace 00112233 with your actual CID)
#ifndef _ADV_PROG_Matrix_06031927_06031927_H_
#define _ADV_PROG_Matrix_06031927_06031927_H_

#include <vector>
#include <iostream>
#include <stdexcept>

namespace adv_prog_cw 
{
    // Step 2: RENAME THE Matrix_06031927 CLASS To Matrix_06031927_00112233
    // (Replace 00112233 with your actual CID)
    template<typename fT>
    class Matrix_06031927 {
    public:
        Matrix_06031927();
        Matrix_06031927(size_t m, size_t n);
        Matrix_06031927(size_t m, size_t n, fT val);
        Matrix_06031927(const Matrix_06031927& M);
        ~Matrix_06031927();

        size_t Rows() const;
        size_t Cols() const;
        void Resize(size_t m, size_t n);
        // accessors M(i,j)
        fT& operator()(size_t m, size_t n);
        const fT& operator()(size_t m, size_t n) const;
        // assignment
        Matrix_06031927& operator=(const Matrix_06031927& M);
        Matrix_06031927& operator=(fT val);

        Matrix_06031927& operator+=(const Matrix_06031927& M);
        Matrix_06031927& operator-=(const Matrix_06031927& M);
        void Identity();
        void Zero();
        void Transposed(Matrix_06031927& RES) const;
        void Out(long digits = 5L) const;

        // Step 3: Additional methods
        Matrix_06031927& operator*=(fT scalar);
        Matrix_06031927& operator/=(fT scalar);
        fT Determinant() const;
        bool Inverse(Matrix_06031927& result) const;
        bool Inverse2(Matrix_06031927& result) const;
        
    private:
        std::vector<std::vector<fT> > data;
        size_t rows, cols;
    
        bool CheckRange(size_t m, size_t n, const char* originator) const;
        bool CheckSizes(const Matrix_06031927& mat, const char* originator) const;
        bool DecomposeLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const;
        bool ParallelLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const;
    };

    // Associated operator declarations
    template<typename fT>
    Matrix_06031927<fT> operator+(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

    template<typename fT>
    Matrix_06031927<fT> operator-(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

    template<typename fT>
    Matrix_06031927<fT> operator*(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

    template<typename fT>
    std::vector<fT> operator*(const Matrix_06031927<fT>& mat, const std::vector<fT>& v);

    // v^T = (v^T M)^T
    template<typename fT>
    Matrix_06031927<fT> operator*(const std::vector<fT>& v, const Matrix_06031927<fT>& M);
} // end namespace

#endif
