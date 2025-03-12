// Step 1: RENAME THE _ADV_PROG_MATRIX_H_ Label To _ADV_PROG_MATRIX_%Your CID Nr%_H_. So, if your CID is 00112233, then the label should be named: _ADV_PROG_MATRIX_00112233_H_
#ifndef _ADV_PROG_MATRIX_H_
#define _ADV_PROG_MATRIX_H_

#include <vector>
#include <iostream>

namespace adv_prog_cw 
{
	// Step 2: RENAME THE MATRIX CLASS To Matrix_%Your CID Nr%. So, if your CID is 00112233, then the class should be named: Matrix_00112233
	template<typename fT>
	class Matrix {
	public:
		Matrix();
		Matrix(size_t m, size_t n);
		Matrix(size_t m, size_t n, fT val);
		Matrix(const Matrix& M);
		~Matrix();

		size_t Rows() const;
		size_t Cols() const;
		void Resize(size_t m, size_t n);
		// accessors M(i,j)
		fT& operator()(size_t m, size_t n);
		const fT& operator()(size_t m, size_t n) const;
		// assignment
		Matrix& operator=(const Matrix& M);
		Matrix& operator=(fT val);

		Matrix& operator+=(const Matrix& M);
		Matrix& operator-=(const Matrix& M);
		void Identity();
		void Zero();
		void Transposed(Matrix& RES) const;
		void Out(long digits = 5L) const;

		// Step 3: Implement the following methods. These will be tested for speed and accuracy for matrices of increasing size and complexity.
		
		// Step 3.1: Multiplication of the matrix by a scalar
		Matrix& operator*=(fT scalar);
		// Step 3.2:  Division of the matrix by a scalar
		Matrix& operator/=(fT scalar);
		// Step 3.3:  A method to compute the determinant of square matrices
		fT Determinant() const;
		// Step 3.4:  A method to compute the inverse of the matrix
		bool Inverse(Matrix& result) const;
		
	private:
		std::vector<std::vector<fT> >  data;
		size_t                   rows, cols;

		bool  CheckRange(size_t m, size_t n, const char* originator) const;
		bool  CheckSizes(const Matrix& mat, const char* originator) const;
	};

	// associated operators
	template<typename fT>
	Matrix<fT>  operator+(const Matrix<fT>& a, const Matrix<fT>& b);

	template<typename fT>
	Matrix<fT>  operator-(const Matrix<fT>& a, const Matrix<fT>& b);

	template<typename fT>
	Matrix<fT>  operator*(const Matrix<fT>& a, const Matrix<fT>& b);

	template<typename fT>
	std::vector<fT>  operator*(const Matrix<fT>& mat, const std::vector<fT>& v);

	// v^T = (v^T M)^T
	template<typename fT>
	Matrix<fT>  operator*(const std::vector<fT>& v, const Matrix<fT>& M);

	// <!-- ==== class Method documentation format ==================== -->

	// constructor (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline Matrix<fT>::Matrix(size_t m, size_t n)
		: data(m, std::vector<fT>(n)), rows(m), cols(n)
	{
	}

	// operator (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline fT& Matrix<fT>::operator()(size_t m, size_t n)
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix<fT>::operator()");
#endif
		return data[m][n];
	}

	// operator (i,j) const
	// ---------------------------------------
	//
	template<typename fT>
	inline const fT& Matrix<fT>::operator()(size_t m, size_t n) const
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix<fT>::operator()");
#endif
		return data[m][n];
	}

	// Rows()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix<fT>::Rows() const { return rows; }

	// Cols()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix<fT>::Cols() const { return cols; }

	// CheckRange(i,j, message )
	// ---------------------------------------
	//
	template<typename fT>
	inline bool Matrix<fT>::CheckRange(size_t m, size_t n,
		const char* originator) const
	{
		if (m >= rows) {
			std::cerr << "\n" << originator << " row index violation, index=" << m;
			std::cerr << " versus, row-max=" << rows << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckRange");
			return false;
		}
		if (n >= cols) {
			std::cerr << "\n" << originator << " column index violation, index=" << n;
			std::cerr << " versus, column-max=" << cols << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckRange");
			return false;
		}
		return true;
	}

	// CheckSizes (i,j,message)
	// ---------------------------------------
	template<typename fT>
	inline bool Matrix<fT>::CheckSizes(const Matrix& mat,
		const char* originator) const
	{
		if (rows != mat.rows) {
			std::cerr << "\n" << originator << " matrices have different sizes; rows1=" << rows;
			std::cerr << " versus, rows2=" << mat.rows << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckSizes");
			return false;
		}
		if (cols != mat.cols) {
			std::cerr << "\n" << originator << " matrices have different sizes; columns1=" << cols;
			std::cerr << " versus, columns2=" << mat.cols << std::endl;
			throw std::length_error("Matrix<double,mn_max>::CheckSizes");
			return false;
		}
		return true;
	}
} // end scope

#endif
