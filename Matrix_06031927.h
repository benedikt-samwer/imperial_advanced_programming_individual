// Step 1: RENAME THE _ADV_PROG_Matrix_06031927_H_ Label To _ADV_PROG_Matrix_06031927_%Your CID Nr%_H_. So, if your CID is 00112233, then the label should be named: _ADV_PROG_Matrix_06031927_00112233_H_
#ifndef _ADV_PROG_Matrix_06031927_06031927_H_
#define _ADV_PROG_Matrix_06031927_06031927_H_

#include <vector>
#include <iostream>

namespace adv_prog_cw 
{
	// Step 2: RENAME THE Matrix_06031927 CLASS To Matrix_06031927_%Your CID Nr%. So, if your CID is 00112233, then the class should be named: Matrix_06031927_00112233
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

		// Step 3: Implement the following methods. These will be tested for speed and accuracy for matrices of increasing size and complexity.
		
		// Step 3.1: Multiplication of the Matrix_06031927 by a scalar
		Matrix_06031927& operator*=(fT scalar);
		// Step 3.2:  Division of the Matrix_06031927 by a scalar
		Matrix_06031927& operator/=(fT scalar);
		// Step 3.3:  A method to compute the determinant of square matrices
		fT Determinant() const;
		// Step 3.4:  A method to compute the inverse of the Matrix_06031927
		bool Inverse(Matrix_06031927& result) const;
		// Step 3.4:  A method to compute the inverse of the Matrix_06031927
		bool Inverse2(Matrix_06031927& result) const;
		
	private:
		std::vector<std::vector<fT> >  data;
		size_t                   rows, cols;
	
		bool  CheckRange(size_t m, size_t n, const char* originator) const;
		bool  CheckSizes(const Matrix_06031927& mat, const char* originator) const;
		bool  DecomposeLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const;
		bool ParallelLU(Matrix_06031927& A, std::vector<size_t>& perm, int& swaps) const;
	};

	// associated operators
	template<typename fT>
	Matrix_06031927<fT>  operator+(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

	template<typename fT>
	Matrix_06031927<fT>  operator-(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

	template<typename fT>
	Matrix_06031927<fT>  operator*(const Matrix_06031927<fT>& a, const Matrix_06031927<fT>& b);

	template<typename fT>
	std::vector<fT>  operator*(const Matrix_06031927<fT>& mat, const std::vector<fT>& v);

	// v^T = (v^T M)^T
	template<typename fT>
	Matrix_06031927<fT>  operator*(const std::vector<fT>& v, const Matrix_06031927<fT>& M);

	// <!-- ==== class Method documentation format ==================== -->

	// constructor (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline Matrix_06031927<fT>::Matrix_06031927(size_t m, size_t n)
		: data(m, std::vector<fT>(n)), rows(m), cols(n)
	{
	}

	// operator (i,j)
	// ---------------------------------------
	//
	template<typename fT>
	inline fT& Matrix_06031927<fT>::operator()(size_t m, size_t n)
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix_06031927<fT>::operator()");
#endif
		return data[m][n];
	}

	// operator (i,j) const
	// ---------------------------------------
	//
	template<typename fT>
	inline const fT& Matrix_06031927<fT>::operator()(size_t m, size_t n) const
	{
#ifndef NDEBUG
		CheckRange(m, n, "Matrix_06031927<fT>::operator()");
#endif
		return data[m][n];
	}

	// Rows()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix_06031927<fT>::Rows() const { return rows; }

	// Cols()
	// ---------------------------------------------
	//
	template<typename fT>
	inline size_t Matrix_06031927<fT>::Cols() const { return cols; }

	// CheckRange(i,j, message )
	// ---------------------------------------
	//
	template<typename fT>
	inline bool Matrix_06031927<fT>::CheckRange(size_t m, size_t n,
		const char* originator) const
	{
		if (m >= rows) {
			std::cerr << "\n" << originator << " row index violation, index=" << m;
			std::cerr << " versus, row-max=" << rows << std::endl;
			throw std::length_error("Matrix_06031927<double,mn_max>::CheckRange");
			return false;
		}
		if (n >= cols) {
			std::cerr << "\n" << originator << " column index violation, index=" << n;
			std::cerr << " versus, column-max=" << cols << std::endl;
			throw std::length_error("Matrix_06031927<double,mn_max>::CheckRange");
			return false;
		}
		return true;
	}

	// CheckSizes (i,j,message)
	// ---------------------------------------
	template<typename fT>
	inline bool Matrix_06031927<fT>::CheckSizes(const Matrix_06031927& mat,
		const char* originator) const
	{
		if (rows != mat.rows) {
			std::cerr << "\n" << originator << " matrices have different sizes; rows1=" << rows;
			std::cerr << " versus, rows2=" << mat.rows << std::endl;
			throw std::length_error("Matrix_06031927<double,mn_max>::CheckSizes");
			return false;
		}
		if (cols != mat.cols) {
			std::cerr << "\n" << originator << " matrices have different sizes; columns1=" << cols;
			std::cerr << " versus, columns2=" << mat.cols << std::endl;
			throw std::length_error("Matrix_06031927<double,mn_max>::CheckSizes");
			return false;
		}
		return true;
	}
} // end scope

#endif
