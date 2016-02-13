/*
 * Matrix.hpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#ifndef MATRIX_HPP_
#define MATRIX_HPP_

class Matrix {

protected:
	double* top_;	// pointer to the matrix
	int m_;			// number of lows of the matrix or lda
	int n_;			// number of columns of the matrix

public:
	// Default constructor
	Matrix();

	// Constructor
	Matrix( const int m, const int n );

	// Copy constructor
	Matrix( const Matrix& T );

	// Destructor
	virtual ~Matrix();

	// Getters
	double* top() { return top_; }
	int m() const { return m_; }
	int n() const { return n_; }

	// Show elements to the standard output
	void Show_all() const;

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Set matrix to the identity matrix
	void Set_Iden();

	// Set matrix to the zero matrix
	void Set_Zero();

	// Assign the value to (i,j) element
	void Set_Val( const int i, const int j, const double val );

	// Operator overload
	Matrix &operator=( const Matrix& T );
	double &operator[]( const int i ) const;
	double &operator()( const int i, const int j ) const;

	// Save matrix elements to the file
	void File_Out( const char* fname );
	void File_Out( const char* fname, const unsigned dig );

};

#endif /* MATRIX_HPP_ */
