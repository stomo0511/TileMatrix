/*
 * BMatrix.hpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#ifndef BMATRIX_HPP_
#define BMATRIX_HPP_

#include "Matrix.hpp"

class BMatrix : public Matrix {

protected:
	int ib_;	// block size

public:
	// Default constructor
	BMatrix();

	// Constructor
	BMatrix( const int m, const int n, const int ib );

	// Copy constructor
	BMatrix( const BMatrix& T );

	// Destructor
	virtual ~BMatrix();

	// Getters
	 int ib() { return ib_; }
};

#endif /* BMATRIX_HPP_ */
