#include "Matrix.h"
#include <cassert>

using namespace std;

typedef int                  int32;
typedef unsigned int         uint32;
typedef double               double64;

namespace adv_prog_cw {
	// default constructor
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>::Matrix()
		: rows(0), cols(0)
	{
	}

	// operator=( DM )
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>& Matrix<fT>::operator=(const Matrix<fT>& mat)
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
	Matrix<fT>& Matrix<fT>::operator=(fT val)
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] = val;

		return *this;
	}

	// copy constructor
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>::Matrix(const Matrix<fT>& mat)
	{
		*this = mat;
	}

	// constructor (i,j, value)
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>::Matrix(size_t m, size_t n, fT val)
		: rows(m), cols(n)
	{
		data.resize(rows, vector<fT>(cols));

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] = val;
	}

	// destructor
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>::~Matrix()
	{
	}

	// operator+=
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>& Matrix<fT>::operator+=(const Matrix<fT>& mat)
	{
		CheckSizes(mat, "Matrix<fT>::operator+=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] += mat.data[i][j];

		return *this;
	}

	// operator-=
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>& Matrix<fT>::operator-=(const Matrix<fT>& mat)
	{
		CheckSizes(mat, "Matrix<fT>::operator-=");
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) data[i][j] -= mat.data[i][j];

		return *this;
	}

	// operator +
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>  operator+(const Matrix<fT>& a, const Matrix<fT>& b)
	{
		Matrix<fT> temp(a);
		temp += b;
		return temp;
	}

	// operator -
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>  operator-(const Matrix<fT>& a, const Matrix<fT>& b)
	{
		Matrix<fT> temp(a);
		temp -= b;
		return temp;
	}

	// operator *
	// ---------------------------------------
	template<typename fT>
	Matrix<fT>  operator*(const Matrix<fT>& a, const Matrix<fT>& b)
	{
#ifndef NDEBUG
		if (a.Cols() != b.Rows()) {
			cout << "\nMatrix<" << typeid(fT).name() << "> operator*: Matrices cannot be multiplied ";
			cout << "because of incompatible sizes (A * B, see matrices below): " << endl;
			a.Out(3L);
			b.Out(3L);
			throw length_error("Matrix<double,mn_max>::operator*");
		}
#endif
		Matrix<fT>  temp(a.Rows(), b.Cols());

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

	// vector^T = Matrix * vector^T
	// ---------------------------------------
	template<typename fT>
	vector<fT>  operator*(const Matrix<fT>& mat, const vector<fT>& vec)
	{
		assert(mat.Cols() == vec.size());

		vector<fT> temp(mat.Rows(), static_cast<fT>(0.0));

		for (size_t i = 0; i < mat.Rows(); i++)
			for (size_t j = 0; j < mat.Cols(); j++) temp[i] += mat(i, j) * vec[j];

		return temp;
	} // end operator*

   // Matrix = vector^T * Matrix
   // ---------------------------------------
	template<typename fT>
	Matrix<fT>  operator*(const vector<fT>& vec, const Matrix<fT>& mat)
	{
		if (vec.size() != mat.Rows()) {
			cerr << "\noperator*: vector cannot be multiplied with matrix";
			cerr << "because of incompatible sizes (v * M): " << endl;
			throw length_error("Matrix<double,mn_max>::operator*");
		}
		Matrix<fT>  temp(vec.size(), mat.Cols());

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
	void Matrix<fT>::Resize(size_t m, size_t n)
	{
		// resize matrix but keep storage as is if the new matrix
		// is smaller than the old matrix
		if (m <= rows && n <= cols) {
			rows = m;
			cols = n;
			return;
		}

		// increase matrix size
		rows = m;
		cols = n;
		data.resize(m);
		for (size_t i = 0; i < m; i++) data[i].resize(n);
	}

	// Identity()
	// ---------------------------------------------
	template<typename fT>
	void Matrix<fT>::Identity()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				if (i == j) data[i][j] = static_cast<fT>(1.0);
				else          data[i][j] = static_cast<fT>(0.0);
	}

	// Zero()
	// ---------------------------------------------
	template<typename fT>
	void Matrix<fT>::Zero()
	{
		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++)
				data[i][j] = static_cast<fT>(0.0);
	}

	// Transposed()
	// ---------------------------------------------
	template<typename fT>
	void Matrix<fT>::Transposed(Matrix& M) const
	{
		M.Resize(cols, rows);

		for (size_t i = 0U; i < rows; i++)
			for (size_t j = 0U; j < cols; j++) M.data[j][i] = data[i][j];
	}

	// Out( digits )
	// ---------------------------------------
	template<typename fT>
	void Matrix<fT>::Out(long digits) const
	{
		std::streamsize  prec;
		cout << "\nMatrix<" << typeid(fT).name() << ">::Out(): m=" << rows << ", n=" << cols << endl;
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

   // ------------------------------------------------------------------------------
   // template instantiations
   // ------------------------------------------------------------------------------

	template class Matrix<double>;

	template Matrix<double>  operator+(
		const Matrix<double>& a,
		const Matrix<double>& b);

	template Matrix<double>  operator-(
		const Matrix<double>& a,
		const Matrix<double>& b);

	template Matrix<double>  operator*(
		const Matrix<double>& a,
		const Matrix<double>& b);
	template
		vector<double>  operator*(const Matrix<double>& mat,
			const vector<double>& vec);

	// extra operators
	template Matrix<double64>
		operator*(const vector<double64>& vec, const Matrix<double64>& mat);
}