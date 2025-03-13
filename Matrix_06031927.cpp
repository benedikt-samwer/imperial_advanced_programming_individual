#include "Matrix_06031927.h"
#include <cassert>

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



	// Scalar Multiplication Operator
    // ---------------------------------------
    // Multiplies each element of the matrix by a scalar.
    // Parameters:
    //   scalar - The scalar value to multiply each element by (type fT)
    // Returns:
    //   Reference to the modified matrix (*this)
    template<typename fT>
    Matrix_06031927<fT>& Matrix_06031927<fT>::operator*=(fT scalar) {
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                data[i][j] *= scalar;
            }
        }
        return *this;
    }

    // Scalar Division Operator
    // ---------------------------------------
    // Divides each element of the matrix by a scalar.
    // Parameters:
    //   scalar - The scalar value to divide each element by (type fT)
    // Returns:
    //   Reference to the modified matrix (*this)
    // Throws:
    //   std::invalid_argument if scalar is zero
    template<typename fT>
    Matrix_06031927<fT>& Matrix_06031927<fT>::operator/=(fT scalar) {
        if (scalar == 0) {
            throw std::invalid_argument("Division by zero in Matrix_06031927::operator/=");
        }
        #pragma omp parallel for collapse(2)
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                data[i][j] /= scalar;
            }
        }
        return *this;
    }

    // LU Decomposition Helper Function
    // ---------------------------------------
    // Performs LU decomposition with partial pivoting on a copy of the matrix.
    // Parameters:
    //   A - Matrix to decompose (modified in-place to store L and U)
    //   perm - Vector storing the row permutation
    //   swaps - Number of row swaps performed (affects determinant sign)
    // Returns:
    //   true if the matrix is non-singular, false if singular
    template<typename fT>
    bool Matrix_06031927<fT>::DecomposeLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const {
        size_t n = A.Rows();
        perm.resize(n);
        for (size_t i = 0; i < n; i++) perm[i] = i;
        swaps = 0;
        for (size_t k = 0; k < n; k++) {
            // Find pivot
            size_t r = k;
            fT max_val = std::abs(A(perm[k], k));
            for (size_t i = k + 1; i < n; i++) {
                fT val = std::abs(A(perm[i], k));
                if (val > max_val) {
                    max_val = val;
                    r = i;
                }
            }
            if (max_val == 0) {
                return false; // Singular matrix
            }
            if (r != k) {
                std::swap(perm[k], perm[r]);
                swaps++;
            }
            fT pivot = A(perm[k], k);
            #pragma omp parallel for
            for (size_t i = k + 1; i < n; i++) {
                fT m = A(perm[i], k) / pivot;
                for (size_t j = k + 1; j < n; j++) {
                    A(perm[i], j) -= m * A(perm[k], j);
                }
                A(perm[i], k) = m;
            }
        }
        return true;
    }

    // Determinant Calculation
    // ---------------------------------------
    // Computes the determinant of a square matrix using LU decomposition.
    // Returns:
    //   The determinant value (type fT); returns 0 if the matrix is singular
    // Throws:
    //   std::invalid_argument if the matrix is not square
    template<typename fT>
    fT Matrix_06031927<fT>::Determinant() const {
        if (rows != cols) {
            throw std::invalid_argument("Matrix must be square for determinant calculation");
        }
        Matrix_06031927 A = *this; // Create a copy since this is const
        std::vector<size_t> perm;
        int swaps;
        if (!DecomposeLU(A, perm, swaps)) {
            return static_cast<fT>(0); // Singular matrix
        }
        fT det = static_cast<fT>(1.0);
        for (size_t i = 0; i < rows; i++) {
            det *= A(perm[i], i);
        }
        if (swaps % 2 == 1) {
            det = -det;
        }
        return det;
    }

    // Matrix Inversion
    // ---------------------------------------
    // Computes the inverse of the matrix and stores it in the result parameter.
    // Parameters:
    //   result - Matrix to store the inverse (resized to match this matrix)
    // Returns:
    //   true if the matrix is invertible, false if not (singular or not square)
    template<typename fT>
    bool Matrix_06031927<fT>::Inverse(Matrix_06031927& result) const {
        if (rows != cols) {
            return false; // Not square
        }
        Matrix_06031927 A = *this; // Create a copy since this is const
        std::vector<size_t> perm;
        int swaps;
        if (!DecomposeLU(A, perm, swaps)) {
            return false; // Singular matrix
        }
        size_t n = rows;
        result.Resize(n, n);
        std::vector<size_t> perm_inv(n);
        for (size_t i = 0; i < n; i++) {
            perm_inv[perm[i]] = i;
        }
        #pragma omp parallel for
        for (size_t j = 0; j < n; j++) {
            // Solve A * x = e_j for column j of the inverse
            std::vector<fT> b(n, static_cast<fT>(0));
            b[perm_inv[j]] = static_cast<fT>(1); // P * e_j
            // Forward substitution: L * y = P * e_j
            std::vector<fT> y(n, static_cast<fT>(0));
            for (size_t i = 0; i < n; i++) {
                fT sum = static_cast<fT>(0);
                for (size_t k = 0; k < i; k++) {
                    sum += A(perm[i], k) * y[k];
                }
                y[i] = b[i] - sum;
            }
            // Backward substitution: U * x = y
            std::vector<fT> x(n, static_cast<fT>(0));
            for (int i = n - 1; i >= 0; i--) {
                fT sum = static_cast<fT>(0);
                for (size_t k = i + 1; k < n; k++) {
                    sum += A(perm[i], k) * x[k];
                }
                x[i] = (y[i] - sum) / A(perm[i], i);
            }
            // Set column j of result
            for (size_t i = 0; i < n; i++) {
                result(i, j) = x[i];
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