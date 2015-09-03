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
	double* top_;		// pointer to the matrix
	unsigned int m_;	// number of lows of the matrix or lda
	unsigned int n_;	// number of columns of the matrix

public:
	// Default constructor
	Matrix();

	// Constructor
	Matrix( const unsigned int m, const unsigned int n );

	// Destructor
	virtual ~Matrix();

	// Getters
	double* top() { return top_; }
	unsigned int m() const { return m_; }
	unsigned int n() const { return n_; }

	// Show elements to the standard output
	void Show_all() const;

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Assign the value to (i,j) element
	void Set_Val( const unsigned int i, const unsigned int j, const double val );

	// Operator overload
	Matrix &operator=( const Matrix& T );
	double &operator[]( const unsigned int i ) const;
	double &operator()( const unsigned int i, const unsigned int j ) const;
};

#endif /* MATRIX_HPP_ */
