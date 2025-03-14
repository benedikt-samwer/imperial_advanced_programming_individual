#include "Matrix_06031927.h"
#include <vector>
#include <cmath>
#include <limits>
#include <algorithm>
#include <stdexcept>
#include <cassert>
#include <cstdlib>
#include <ctime>
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



	// Scalar Multiplication Operator
    // ---------------------------------------
    // Multiplies each element of the matrix by a scalar.
    // Parameters:
    //   scalar - The scalar value to multiply each element by (type fT)
    // Returns:
    //   Reference to the modified matrix (*this)
	template<typename fT>
	Matrix_06031927<fT>& Matrix_06031927<fT>::operator*=(fT scalar) {
		if (rows * cols > 10000) {  // Threshold: 100x100 matrix
			#pragma omp parallel for collapse(2)
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					data[i][j] *= scalar;
				}
			}
		} else {
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					data[i][j] *= scalar;
				}
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
		if (rows * cols > 10000) {  // Threshold: 100x100 matrix
			#pragma omp parallel for collapse(2)
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					data[i][j] /= scalar;
				}
			}
		} else {
			for (size_t i = 0; i < rows; i++) {
				for (size_t j = 0; j < cols; j++) {
					data[i][j] /= scalar;
				}
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
	bool Matrix_06031927<fT>::ParallelLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps, bool& earlyExit, double& log_det, int& sign_det) const 
	{
		size_t n = A.Rows();
		perm.resize(n);
		for (size_t i = 0; i < n; i++) {
			perm[i] = i;
		}
		swaps = 0;
		earlyExit = false;
		log_det = 0.0;
		sign_det = 1; // Start with positive sign.

		const size_t blockSize = 128;      // Tuned for cache efficiency
		const double INF_THRESHOLD = 710.0;  // ln(max_double) ~709.78

		for (size_t k = 0; k < n; k++) {
			// --- Pivot Selection ---
			size_t pivotRow = k;
			fT max_val = std::abs(A(perm[k], k));
			#pragma omp parallel
			{
				size_t local_pivotRow = k;
				fT local_max_val = max_val;
				#pragma omp for nowait
				for (size_t i = k + 1; i < n; i++) {
					fT cur_val = std::abs(A(perm[i], k));
					if (cur_val > local_max_val) {
						local_max_val = cur_val;
						local_pivotRow = i;
					}
				}
				#pragma omp critical
				{
					if (local_max_val > max_val) {
						max_val = local_max_val;
						pivotRow = local_pivotRow;
					}
				}
			}
			if (max_val == static_cast<fT>(0))
				return false;  // Singular matrix

			if (std::isinf(max_val)) {
				earlyExit = true;
				log_det += std::log(max_val);
				break;
			}

			// --- Row Swap ---
			if (pivotRow != k) {
				std::swap(perm[k], perm[pivotRow]);
				swaps++;
				sign_det = -sign_det;  // Row swap flips the sign
			}

			// --- Process the Pivot ---
			fT pivot = A(perm[k], k);
			if (pivot < 0)
				sign_det = -sign_det;  // Pivot sign affects determinant

			log_det += std::log(std::abs(pivot));
			if (log_det > INF_THRESHOLD) {
				earlyExit = true;
				break;
			}

			// --- LU Update ---
			#pragma omp parallel for schedule(static)
			for (size_t i = k + 1; i < n; i++) {
				fT m = A(perm[i], k) / pivot;
				A(perm[i], k) = m;
				for (size_t j = k + 1; j < n; j += blockSize) {
					size_t j_end = std::min(j + blockSize, n);
					#pragma omp simd
					for (size_t jj = j; jj < j_end; jj++) {
						A(perm[i], jj) -= m * A(perm[k], jj);
					}
				}
			}
		}
		return true;
	}

	template<typename fT>
	fT Matrix_06031927<fT>::Determinant() const {
		if (rows != cols) {
			throw std::invalid_argument("Matrix must be square for determinant calculation");
		}

		Matrix_06031927 A = *this;  // Work on a copy (LU is in-place)
		std::vector<size_t> perm;
		int swaps = 0;
		bool earlyExit = false;
		double log_det = 0.0;
		int sign_det = 1;

		if (!ParallelLU(A, perm, swaps, earlyExit, log_det, sign_det))
			return static_cast<fT>(0);

		// Handle overflow case
		if (earlyExit || log_det > 710.0) {
			return (sign_det > 0) ? std::numeric_limits<fT>::infinity() 
								: -std::numeric_limits<fT>::infinity();
		}

		// Compute determinant as sign_det * exp(log_det)
		return static_cast<fT>(sign_det) * std::exp(log_det);
	}


	template<typename fT>
	bool Matrix_06031927<fT>::ParallelLUforinverse(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const {
		size_t n = A.Rows();
		perm.resize(n);
		for (size_t i = 0; i < n; i++) {
			perm[i] = i;
		}
		swaps = 0;

		for (size_t k = 0; k < n; k++) {
			// --- Parallel Pivot Selection ---
			size_t pivotRow = k;
			fT max_val = std::abs(A(perm[k], k));
			#pragma omp parallel
			{
				size_t local_pivotRow = k;
				fT local_max_val = max_val;
				#pragma omp for nowait
				for (size_t i = k + 1; i < n; i++) {
					fT cur_val = std::abs(A(perm[i], k));
					if (cur_val > local_max_val) {
						local_max_val = cur_val;
						local_pivotRow = i;
					}
				}
				#pragma omp critical
				{
					if (local_max_val > max_val) {
						max_val = local_max_val;
						pivotRow = local_pivotRow;
					}
				}
			}
			if (max_val == static_cast<fT>(0)) {
				return false;  // Singular matrix
			}

			// --- Swap pivot row ---
			if (pivotRow != k) {
				std::swap(perm[k], perm[pivotRow]);
				swaps++;
			}
			fT pivot = A(perm[k], k);

			// --- Update Trailing Submatrix ---
			#pragma omp parallel for schedule(static)
			for (size_t i = k + 1; i < n; i++) {
				fT m = A(perm[i], k) / pivot;
				A(perm[i], k) = m;
				#pragma omp simd
				for (size_t j = k + 1; j < n; j++) {
					A(perm[i], j) -= m * A(perm[k], j);
				}
			}
		}
		return true;
	}

	template<typename fT>
	bool Matrix_06031927<fT>::Inverse(Matrix_06031927& result) const {
		if (rows != cols) return false;  // Must be square
		size_t n = rows;

		// Step 1: Make a copy for in-place LU decomposition.
		Matrix_06031927 A = *this;
		std::vector<size_t> perm(n);
		for (size_t i = 0; i < n; i++) {
			perm[i] = i;
		}
		int swaps;
		// Use the fully computed, parallel LU for inversion.
		if (!ParallelLUforinverse(A, perm, swaps))
			return false; // Singular matrix

		// Step 2: Compute the inverse permutation vector.
		std::vector<size_t> perm_inv(n);
		for (size_t i = 0; i < n; i++) {
			perm_inv[perm[i]] = i;
		}

		// Step 3: Resize the result matrix.
		result.Resize(n, n);

		// Step 4: Compute inverse columns in parallel.
		#pragma omp parallel for schedule(dynamic)
		for (size_t j = 0; j < n; j++) {
			// Allocate per-thread local vectors.
			// (Allocation here ensures they are private to each column solve.)
			std::vector<fT> b(n, 0.0);
			std::vector<fT> y(n, 0.0);
			std::vector<fT> x(n, 0.0);

			// Compute the right-hand side: b = P * e_j.
			b[perm_inv[j]] = 1.0;

			// Forward substitution: Solve L * y = b.
			for (size_t i = 0; i < n; i++) {
				fT sum = 0.0;
				size_t row_index = perm[i];
				#pragma omp simd reduction(+:sum)
				for (size_t k = 0; k < i; k++) {
					sum += A(row_index, k) * y[k];
				}
				y[i] = b[i] - sum;
			}

			// Backward substitution: Solve U * x = y.
			for (int i = static_cast<int>(n) - 1; i >= 0; i--) {
				fT sum = 0.0;
				size_t row_index = perm[i];
				#pragma omp simd reduction(+:sum)
				for (size_t k = i + 1; k < n; k++) {
					sum += A(row_index, k) * x[k];
				}
				x[i] = (y[i] - sum) / A(row_index, i);
			}

			// Store the computed column in the result.
			// If your result matrix is stored contiguously, consider using std::copy.
			for (size_t i = 0; i < n; i++) {
				result(i, j) = x[i];
			}
		}
		return true;
	}

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