/*
 * BMatrix.cpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include "BMatrix.hpp"

using namespace std;

/**
 * Default constructor
 */
BMatrix::BMatrix() : Matrix()
{
	#ifdef DEBUG
	cout << "BMatrix()\n";
	#endif

	ib_ = 0;
}

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
BMatrix::BMatrix( const int m, const int n, const int ib ) : Matrix(m,n)
{
	#ifdef DEBUG
	cout << "BMatrix(m,n,ib)\n";
	#endif

	assert( ib > 0 );
	assert( ib < m );

	ib_ = ib;
}

/**
 * Destructor
 *
 */
BMatrix::~BMatrix()
{
	#ifdef DEBUG
	cout << "~BMatrix()\n";
	#endif
}
