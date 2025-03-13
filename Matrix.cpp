#include "Matrix_06031927.h"
#include <cassert>
#include <vector>
#include <cmath>
#include <omp.h>
#include <stdexcept>
#include <cmath>
#include <cstdlib>
#include <vector>
#include <stdexcept>
#include <limits>
#include <cassert>
#include <ctime>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>


using namespace std;

typedef int                  int32;
typedef unsigned int         uint32;
typedef double               double64;

namespace adv_prog_cw {
	// default constructor
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>::Matrix_06031927()
		: rows(0), cols(0)
	{
	}

	// operator=( DM )
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>& Matrix_06031927<fT>::operator=(const Matrix_06031927<fT>& mat)
	{
		if (&mat != this) {
			rows = mat.rows;
			cols = mat.cols;

			data.resize(rows, vector<fT>(cols));

			for (size_t i = 0U; i < rows; i++)
				for (size_t j = 0U; j < cols; j++)
					data[i][j] = mat.data[i][j];
		}

		return *this;
	}

	template<typename fT>
	Matrix_06031927<fT>& Matrix_06031927<fT>::operator=(fT val)
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] = val;

		return *this;
	}

	// copy constructor
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>::Matrix_06031927(const Matrix_06031927<fT>& mat)
	{
		*this = mat;
	}

	// constructor (i,j, value)
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>::Matrix_06031927(size_t m, size_t n, fT val)
		: rows(m), cols(n)
	{
		data.resize(rows, vector<fT>(cols));

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] = val;
	}

	// destructor
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>::~Matrix_06031927()
	{
	}

	// operator+=
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>& Matrix_06031927<fT>::operator+=(const Matrix_06031927<fT>& mat)
	{
		CheckSizes(mat, "Matrix_06031927<fT>::operator+=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] += mat.data[i][j];

		return *this;
	}

	// operator-=
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>& Matrix_06031927<fT>::operator-=(const Matrix_06031927<fT>& mat)
	{
		CheckSizes(mat, "Matrix_06031927<fT>::operator-=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] -= mat.data[i][j];

		return *this;
	}

	// operator +
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>  operator+(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b)
	{
		Matrix_06031927<fT> temp(a);
		temp += b;
		return temp;
	}

	// operator -
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>  operator-(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b)
	{
		Matrix_06031927<fT> temp(a);
		temp -= b;
		return temp;
	}

	// operator *
	// ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>  operator*(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b)
	{
#ifndef NDEBUG
		if (a.Cols() != b.Rows()) {
			cout << "\nMatrix_06031927<" << typeid(fT).name() << "> operator*: Matrices cannot be multiplied ";
			cout << "because of incompatible sizes (A * B, see matrices below): " << endl;
			a.Out(3L);
			b.Out(3L);
			throw length_error("Matrix_06031927<double,mn_max>::operator*");
		}
#endif
		Matrix_06031927<fT>  temp(a.Rows(), b.Cols());

		for (size_t i = 0U; i < a.Rows(); i++)
			for (size_t j = 0U; j < b.Cols(); j++) {
				temp(i, j) = static_cast<fT>(0.0);
				for (size_t k = 0U; k < b.Rows(); k++)
					temp(i, j) += a(i, k) * b(k, j);
			}

		return temp;
	}

	// -------------------------------------------------------------------------------------------
	// OPERATOR FUNCTIONS
	// -------------------------------------------------------------------------------------------

	// vector^T = Matrix_06031927 * vector^T
	// ---------------------------------------
	template<typename fT>
	vector<fT>  operator*(const Matrix_06031927<fT>& mat, const vector<fT>& vec)
	{
		assert(mat.Cols() == vec.size());

		vector<fT> temp(mat.Rows(), static_cast<fT>(0.0));

		for (size_t i = 0; i < mat.Rows(); i++)
			for (size_t j = 0; j < mat.Cols(); j++) temp[i] += mat(i, j) * vec[j];

		return temp;
	} // end operator*

   // Matrix_06031927 = vector^T * Matrix_06031927
   // ---------------------------------------
	template<typename fT>
	Matrix_06031927<fT>  operator*(const vector<fT>& vec, const Matrix_06031927<fT>& mat)
	{
		if (vec.size() != mat.Rows()) {
			cerr << "\noperator*: vector cannot be multiplied with Matrix_06031927";
			cerr << "because of incompatible sizes (v * M): " << endl;
			throw length_error("Matrix_06031927<double,mn_max>::operator*");
		}
		Matrix_06031927<fT>  temp(vec.size(), mat.Cols());

		for (size_t i = 0U; i < vec.size(); i++)
			for (size_t j = 0U; j < mat.Cols(); j++) {
				temp(i, j) = static_cast<fT>(0.0);
				for (size_t k = 0U; k < mat.Rows(); k++)
					temp(i, j) += vec[k] * mat(k, j);
			}

		return temp;
	} // end operator

   // Resize()
   // ---------------------------------------------
	template<typename fT>
	void Matrix_06031927<fT>::Resize(size_t m, size_t n)
	{
		// resize Matrix_06031927 but keep storage as is if the new Matrix_06031927
		// is smaller than the old Matrix_06031927
		if (m <= rows && n <= cols) {
			rows = m;
			cols = n;
			return;
		}

		// increase Matrix_06031927 size
		rows = m;
		cols = n;
		data.resize(m);
		for (size_t i = 0; i < m; i++) data[i].resize(n);
	}

	// Identity()
	// ---------------------------------------------
	template<typename fT>
	void Matrix_06031927<fT>::Identity()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				if (i == j) data[i][j] = static_cast<fT>(1.0);
				else          data[i][j] = static_cast<fT>(0.0);
	}

	// Zero()
	// ---------------------------------------------
	template<typename fT>
	void Matrix_06031927<fT>::Zero()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] = static_cast<fT>(0.0);
	}

	// Transposed()
	// ---------------------------------------------
	template<typename fT>
	void Matrix_06031927<fT>::Transposed(Matrix_06031927& M) const
	{
		M.Resize(cols, rows);

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) M.data[j][i] = data[i][j];
	}

	// Out( digits )
	// ---------------------------------------
	template<typename fT>
	void Matrix_06031927<fT>::Out(long digits) const
	{
		std::streamsize  prec;
		cout << "\nMatrix_06031927<" << typeid(fT).name() << ">::Out(): m=" << rows << ", n=" << cols << endl;
		if (digits != 0U) {
			cout.setf(ios::scientific);
			prec = cout.precision(digits);
		}
		size_t row_break, split_after(10U);

		for (size_t i = 0; i < rows; i++)
		{
			row_break = 1;
			for (size_t j = 0; j < cols; j++, row_break++)
			{
				if (data[i][j] >= 0.) cout << " ";
				cout << data[i][j] << " ";
				if (row_break == split_after)
				{
					cout << endl;
					row_break = 0U;
				}
			}
			cout << endl;
		}

		if (digits != 0U) {
			cout.unsetf(ios::scientific);
			cout.precision(prec);
		}

		cout << endl;
	} // end Out()


    // ------------------------------------------------------------------
    // Step 3.1: Multiplication of the Matrix by a Scalar
    // ------------------------------------------------------------------
    // Multiplies each element of the matrix by a scalar.
    // Parameters:
    //   - scalar: the value to multiply each element by.
    // Returns:
    //   - Reference to the modified matrix.
    template<typename fT>
    Matrix_06031927<fT>& Matrix_06031927<fT>::operator*=(fT scalar) {
        // Parallelize the outer loop if OpenMP is enabled.
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(scalar)
#endif
        for (size_t i = 0U; i < rows; i++) {
            for (size_t j = 0U; j < cols; j++) {
                data[i][j] *= scalar;
            }
        }
        return *this;
    }

    // ------------------------------------------------------------------
    // Step 3.2: Division of the Matrix by a Scalar
    // ------------------------------------------------------------------
    // Divides each element of the matrix by a scalar.
    // Parameters:
    //   - scalar: the value to divide each element by.
    // Returns:
    //   - Reference to the modified matrix.
    // Throws:
    //   - std::runtime_error if scalar equals zero.
    template<typename fT>
    Matrix_06031927<fT>& Matrix_06031927<fT>::operator/=(fT scalar) {
        if (scalar == static_cast<fT>(0)) {
            throw std::runtime_error("Division by zero in Matrix_06031927::operator/=");
        }
#ifdef _OPENMP
        #pragma omp parallel for default(none) shared(scalar)
#endif
        for (size_t i = 0U; i < rows; i++) {
            for (size_t j = 0U; j < cols; j++) {
                data[i][j] /= scalar;
            }
        }
        return *this;
    }

    // ------------------------------------------------------------------
    // Step 3.3: Compute the Determinant of a Square Matrix
    // ------------------------------------------------------------------
    // Computes and returns the determinant of a square matrix.
    // Returns:
    //   - Determinant value of the matrix.
    // Throws:
    //   - std::length_error if the matrix is not square.

    //---------------------- Helper Functions ----------------------//

    // ExtractBlock: extracts a submatrix from A starting at (rowStart, colStart)
    // with dimensions (numRows x numCols).
    template<typename fT>
    Matrix_06031927<fT> ExtractBlock(const Matrix_06031927<fT>& A, size_t rowStart, size_t colStart, size_t numRows, size_t numCols) {
        Matrix_06031927<fT> block(numRows, numCols, 0);
        for (size_t i = 0; i < numRows; i++) {
            for (size_t j = 0; j < numCols; j++) {
                block(i, j) = A(rowStart + i, colStart + j);
            }
        }
        return block;
    }

    // SetBlock: sets the contents of B into matrix A starting at (rowStart, colStart).
    template<typename fT>
    void SetBlock(Matrix_06031927<fT>& A, const Matrix_06031927<fT>& B, size_t rowStart, size_t colStart) {
        size_t numRows = B.Rows();
        size_t numCols = B.Cols();
        for (size_t i = 0; i < numRows; i++) {
            for (size_t j = 0; j < numCols; j++) {
                A(rowStart + i, colStart + j) = B(i, j);
            }
        }
    }

    // ParallelMultiply: performs matrix multiplication C = A * B using OpenMP.
    template<typename fT>
    void ParallelMultiply(const Matrix_06031927<fT>& A, const Matrix_06031927<fT>& B, Matrix_06031927<fT>& C) {
        size_t m = A.Rows();
        size_t p = A.Cols();
        size_t n = B.Cols();
        C.Resize(m, n);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(m); i++) {
            for (size_t j = 0; j < n; j++) {
                fT sum = 0;
                for (size_t k = 0; k < p; k++) {
                    sum += A(i, k) * B(k, j);
                }
                C(i, j) = sum;
            }
        }
    }

    // MatrixAdd: returns the matrix sum C = A + B.
    template<typename fT>
    Matrix_06031927<fT> MatrixAdd(const Matrix_06031927<fT>& A, const Matrix_06031927<fT>& B) {
        size_t m = A.Rows(), n = A.Cols();
        Matrix_06031927<fT> C(m, n, 0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(m); i++) {
            for (size_t j = 0; j < n; j++) {
                C(i, j) = A(i, j) + B(i, j);
            }
        }
        return C;
    }

    // MatrixSubtract: returns the matrix difference C = A - B.
    template<typename fT>
    Matrix_06031927<fT> MatrixSubtract(const Matrix_06031927<fT>& A, const Matrix_06031927<fT>& B) {
        size_t m = A.Rows(), n = A.Cols();
        Matrix_06031927<fT> C(m, n, 0);
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(m); i++) {
            for (size_t j = 0; j < n; j++) {
                C(i, j) = A(i, j) - B(i, j);
            }
        }
        return C;
    }

    //------------------- Recursive Block Determinant -------------------//

    // BlockDeterminant recursively computes the determinant of a square matrix A.
    // Parameters:
    //   A   - The matrix whose determinant is to be computed.
    //   det - A reference where the computed determinant will be stored.
    // Returns:
    //   true if the determinant is computed successfully, false otherwise.
    template<typename fT>
    bool BlockDeterminant(const Matrix_06031927<fT>& A, fT &det) {
        size_t n = A.Rows();
        const fT epsilon = static_cast<fT>(1e-12);
        
        // Base case: 1x1 matrix.
        if(n == 1) {
            det = A(0, 0);
            return true;
        }
        // Base case: 2x2 matrix using the formula:
        // det = a11*a22 - a12*a21.
        if(n == 2) {
            det = A(0,0) * A(1,1) - A(0,1) * A(1,0);
            return true;
        }
        
        // Partition A into four blocks.
        size_t k = n / 2;       // Size of block A11.
        size_t r = n - k;       // Size of the remaining blocks.
        Matrix_06031927<fT> A11 = ExtractBlock(A, 0, 0, k, k);
        Matrix_06031927<fT> A12 = ExtractBlock(A, 0, k, k, r);
        Matrix_06031927<fT> A21 = ExtractBlock(A, k, 0, r, k);
        Matrix_06031927<fT> A22 = ExtractBlock(A, k, k, r, r);
        
        // Compute determinant of A11 recursively.
        fT detA11;
        if (!BlockDeterminant(A11, detA11))
            return false;
        // If A11 is nearly singular, the determinant is (practically) zero.
        if (std::abs(detA11) < epsilon) {
            det = 0;
            return true;
        }
        
        // To form the Schur complement, we need A11^{-1}.
        // We use our previously defined BlockInverse method.
        Matrix_06031927<fT> A11_inv;
        if (!BlockInverse(A11, A11_inv)) {  // If A11 is singular.
            det = 0;
            return true;
        }
        
        // Compute the product A21 * A11_inv * A12.
        Matrix_06031927<fT> temp;
        ParallelMultiply(A11_inv, A12, temp);
        Matrix_06031927<fT> prod;
        ParallelMultiply(A21, temp, prod);
        
        // Compute the Schur complement: S = A22 - A21*A11^{-1}A12.
        Matrix_06031927<fT> S = MatrixSubtract(A22, prod);
        
        // Recursively compute the determinant of the Schur complement.
        fT detS;
        if (!BlockDeterminant(S, detS))
            return false;
        
        // The determinant of A is given by:
        // det(A) = det(A11) * det(S)
        det = detA11 * detS;
        return true;
    }

    //------------------- Public Determinant() Method -------------------//

    template<typename fT>
    bool Matrix_06031927<fT>::Determinant(fT& det) const {
        // Determinant is only defined for square matrices.
        if (rows != cols)
            return false;
        return BlockDeterminant(*this, det);
    }

    //------------------- Sample Implementations for Other Methods -------------------//

    template<typename fT>
    Matrix_06031927<fT>::Matrix_06031927(size_t m, size_t n, fT val)
        : rows(m), cols(n), data(m, std::vector<fT>(n, val))
    {
    }

    template<typename fT>
    size_t Matrix_06031927<fT>::Rows() const { return rows; }

    template<typename fT>
    size_t Matrix_06031927<fT>::Cols() const { return cols; }

    template<typename fT>
    void Matrix_06031927<fT>::Resize(size_t m, size_t n) {
        rows = m;
        cols = n;
        data.resize(m);
        for (size_t i = 0; i < m; i++)
            data[i].resize(n);
    }

    template<typename fT>
    fT& Matrix_06031927<fT>::operator()(size_t i, size_t j) {
        return data[i][j];
    }

    template<typename fT>
    const fT& Matrix_06031927<fT>::operator()(size_t i, size_t j) const {
        return data[i][j];
    }

    // Destructor and copy constructor would be implemented as needed.

    // ------------------------------------------------------------------
    // Step 3.4: Compute the Inverse of the Matrix
    // ------------------------------------------------------------------
    // Computes the inverse of the matrix and stores the result in 'result'.
    // Parameters:
    //   - result: a Matrix_06031927 object where the inverse will be stored.
    // Returns:
    //   - true if the matrix is invertible, false otherwise.
    // Note:
    //   - This method uses the Gauss-Jordan elimination with partial pivoting.

    
    // -------------------- Recursive Block Inversion -------------------- //

    // This function attempts to compute the inverse of matrix A (which must be square)
    // using block inversion with the Schur complement.
    // Returns true on success, false if a singular sub-block is encountered.
    template<typename fT>
    bool BlockInverse(const Matrix_06031927<fT>& A, Matrix_06031927<fT>& A_inv) {
        size_t n = A.Rows();
        // Base case: 1x1 matrix.
        if(n == 1) {
            const fT eps = static_cast<fT>(1e-12);
            if (std::abs(A(0,0)) < eps)
                return false;
            A_inv.Resize(1,1);
            A_inv(0,0) = static_cast<fT>(1) / A(0,0);
            return true;
        }

        // Partition A into 4 blocks.
        size_t k = n / 2;         // size of A11 (rows and columns)
        size_t r = n - k;         // remaining rows/columns for A22

        Matrix_06031927<fT> A11 = ExtractBlock(A, 0, 0, k, k);
        Matrix_06031927<fT> A12 = ExtractBlock(A, 0, k, k, r);
        Matrix_06031927<fT> A21 = ExtractBlock(A, k, 0, r, k);
        Matrix_06031927<fT> A22 = ExtractBlock(A, k, k, r, r);

        // Invert A11.
        Matrix_06031927<fT> A11_inv;
        if (!BlockInverse(A11, A11_inv))
            return false;  // A11 singular

        // Compute intermediate: M = A21 * A11_inv.
        Matrix_06031927<fT> M;
        ParallelMultiply(A21, A11_inv, M);

        // Compute Schur complement: S = A22 - M * A12.
        Matrix_06031927<fT> M_A12;
        ParallelMultiply(M, A12, M_A12);
        Matrix_06031927<fT> S = MatrixSubtract(A22, M_A12);

        // Invert S.
        Matrix_06031927<fT> S_inv;
        if (!BlockInverse(S, S_inv))
            return false;  // Schur complement singular

        // Now compute the blocks of A_inv.
        // B11 = A11_inv + A11_inv * A12 * S_inv * A21 * A11_inv.
        Matrix_06031927<fT> A12_Sinv;
        ParallelMultiply(A12, S_inv, A12_Sinv);
        Matrix_06031927<fT> temp1;
        ParallelMultiply(A11_inv, A12_Sinv, temp1);  // temp1 = A11_inv * A12 * S_inv
        Matrix_06031927<fT> A21_A11inv;
        ParallelMultiply(A21, A11_inv, A21_A11inv);
        Matrix_06031927<fT> temp2;
        ParallelMultiply(temp1, A21, temp2); // temp2 = A11_inv * A12 * S_inv * A21
        Matrix_06031927<fT> B11;
        {
            Matrix_06031927<fT> correction;
            ParallelMultiply(temp2, A11_inv, correction); // correction = A11_inv * A12 * S_inv * A21 * A11_inv
            B11 = MatrixAdd(A11_inv, correction);
        }

        // B12 = -A11_inv * A12 * S_inv.
        Matrix_06031927<fT> B12;
        {
            Matrix_06031927<fT> temp;
            ParallelMultiply(A11_inv, A12_Sinv, temp);  // temp = A11_inv * A12 * S_inv
            B12 = temp;
            size_t mB = B12.Rows(), nB = B12.Cols();
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < static_cast<int>(mB); i++) {
                for (size_t j = 0; j < nB; j++) {
                    B12(i,j) = -B12(i,j);
                }
            }
        }

        // B21 = -S_inv * A21 * A11_inv.
        Matrix_06031927<fT> B21;
        {
            Matrix_06031927<fT> temp;
            ParallelMultiply(A21, A11_inv, temp);  // temp = A21 * A11_inv
            ParallelMultiply(S_inv, temp, B21);      // B21 = S_inv * (A21 * A11_inv)
            size_t mB = B21.Rows(), nB = B21.Cols();
            #pragma omp parallel for schedule(static)
            for (int i = 0; i < static_cast<int>(mB); i++) {
                for (size_t j = 0; j < nB; j++) {
                    B21(i,j) = -B21(i,j);
                }
            }
        }

        // B22 = S_inv.
        Matrix_06031927<fT> B22 = S_inv; 

        // Combine blocks into A_inv.
        A_inv.Resize(n, n);
        // Top-left block
        SetBlock(A_inv, B11, 0, 0);
        // Top-right block
        SetBlock(A_inv, B12, 0, k);
        // Bottom-left block
        SetBlock(A_inv, B21, k, 0);
        // Bottom-right block
        SetBlock(A_inv, B22, k, k);

        return true;
    }

    // -------------------- Inverse() Method Implementation -------------------- //

    template<typename fT>
    bool Matrix_06031927<fT>::Inverse(Matrix_06031927<fT>& result) const {
        // The matrix must be square.
        if (rows != cols)
            return false;
        return BlockInverse(*this, result);
    }
    template<typename fT>
    bool Matrix_06031927<fT>::Inverse2(Matrix_06031927<fT>& result) const {
        // The matrix must be square for an inverse to exist.
        if (rows != cols)
            return false;
        
        const size_t n = rows; // dimension of the square matrix
        // Create an augmented matrix of size n x 2n.
        // Left half: a copy of the original matrix.
        // Right half: the identity matrix.
        std::vector<std::vector<fT>> aug(n, std::vector<fT>(2 * n, static_cast<fT>(0)));
        for (size_t i = 0; i < n; i++) {
            for (size_t j = 0; j < n; j++) {
                aug[i][j] = data[i][j];
            }
            // Set the identity matrix portion.
            aug[i][n + i] = static_cast<fT>(1);
        }
        
        // Define an epsilon threshold for pivot comparisons.
        const fT epsilon = static_cast<fT>(1e-12);
        
        // Perform Gauss-Jordan elimination to transform [A | I] into [I | A^-1].
        for (size_t i = 0; i < n; i++) {
            // Partial pivoting: find the row with the maximum absolute value in the current column.
            size_t pivotRow = i;
            fT maxVal = std::abs(aug[i][i]);
            for (size_t r = i + 1; r < n; r++) {
                fT val = std::abs(aug[r][i]);
                if (val > maxVal) {
                    maxVal = val;
                    pivotRow = r;
                }
            }
            // If the pivot is too small, the matrix is singular.
            if (maxVal < epsilon)
                return false;
            
            // Swap the current row with the pivot row if necessary.
            if (pivotRow != i)
                std::swap(aug[i], aug[pivotRow]);
            
            // Normalize the pivot row so that the pivot element becomes 1.
            fT pivot = aug[i][i];
            // Parallelize the division of the entire pivot row.
            #pragma omp parallel for schedule(static)
            for (int k = static_cast<int>(i); k < static_cast<int>(2 * n); k++) {
                aug[i][k] /= pivot;
            }
            
            // Eliminate the current pivot column elements in all other rows.
            #pragma omp parallel for schedule(static)
            for (int j = 0; j < static_cast<int>(n); j++) {
                if (j == static_cast<int>(i))
                    continue;
                fT factor = aug[j][i];
                for (int k = static_cast<int>(i); k < static_cast<int>(2 * n); k++) {
                    aug[j][k] -= factor * aug[i][k];
                }
            }
        }
        
        // At this stage, the left half of the augmented matrix should be the identity matrix,
        // and the right half is the computed inverse.
        result.Resize(n, n);
        // Copy the inverse (right half) into the result matrix.
        #pragma omp parallel for schedule(static)
        for (int i = 0; i < static_cast<int>(n); i++) {
            for (int j = 0; j < static_cast<int>(n); j++) {
                result(i, j) = aug[i][j + n];
            }
        }
        
        return true;
    }

    // Existing implementations continue here (operator+, operator-, etc.)...

    // Template instantiations (updated to include new methods implicitly)
    template class Matrix_06031927<double>;

    template Matrix_06031927<double>  operator+(
        const Matrix_06031927<double>& a,
        const Matrix_06031927<double>& b);

    template Matrix_06031927<double>  operator-(
        const Matrix_06031927<double>& a,
        const Matrix_06031927<double>& b);

    template Matrix_06031927<double>  operator*(
        const Matrix_06031927<double>& a,
        const Matrix_06031927<double>& b);

    template
        vector<double>  operator*(const Matrix_06031927<double>& mat,
            const vector<double>& vec);

    template Matrix_06031927<double64>
        operator*(const vector<double64>& vec, const Matrix_06031927<double64>& mat);

}