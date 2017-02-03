/*
 * Matrix.cpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cassert>
#include "Matrix.hpp"

using namespace std;

/**
 * Default constructor
 */
Matrix::Matrix()
{
	#ifdef DEBUG
	cout << "Matrix()\n";
	#endif

	m_ = n_ = 0;
	top_ = NULL;
}

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
Matrix::Matrix( const int m, const int n )
{
	#ifdef DEBUG
	cout << "Matrix(m,n)\n";
	#endif

	assert( m > 0 );
	assert( n > 0 );

	m_ = m;
	n_ = n;

	try
	{
		top_ = new double[ m_ * n_ ];
	}
	catch (char *eb)
	{
		cerr << "Can't allocate memory space for Matrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}
}

/**
 * Copy Constructor
 *
 * @param T Matrix object
 */
Matrix::Matrix( const Matrix& T )
{
	#ifdef DEBUG
	cout << "Matrix(m,n)\n";
	#endif

	assert( T.top_ != NULL );

	m_ = T.m_;
	n_ = T.n_;

	try
	{
		top_ = new double[ m_ * n_ ];
	}
	catch (char *eb)
	{
		cerr << "Can't allocate memory space for Matrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}

	for (int i=0; i< m_ * n_; i++)
		top_[i] = T.top_[i];
}

/**
 * Destructor
 *
 */
Matrix::~Matrix()
{
	#ifdef DEBUG
	cout << "~Matrix()\n";
	#endif

	delete [] top_;
}

/**
 * Show elements to the standard output
 */
void Matrix::Show_all() const
{
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			cout << top_[ i + j * m_ ] << ", ";
		}
		cout << endl;
	}
	cout << endl;
}

/**
 * Assign random numbers to the elements
 *
 * @param seed Seed of random number generator
 */
void Matrix::Set_Rnd( const unsigned seed )
{
	assert( seed >= 0 );

	srand(seed);
	for (int i = 0; i < m_ * n_; i++)
		top_[i] = (double)rand() / RAND_MAX;
}

/**
 * Set matrix to the identity matrix
 */
void Matrix::Set_Iden()
{
	for (int i=0; i<m_; i++)
		for (int j=0; j<n_; j++)
			top_[ i + j*m_ ] = ( i == j ) ? (double)(1) : (double)(0);
}

/**
 * Set matrix to the zero matrix
 */
void Matrix::Set_Zero()
{
	for (int i=0; i< m_ * n_; i++)
		top_[i] = (double)0;
}

/**
 * Assign the value to (i,j) element
 *
 * @param i vertical index of the element
 * @param j horizontal index of the element
 * @param val element value
 */
void Matrix::Set_Val( const int i, const int j, const double val )
{
	assert( i >= 0 );	assert( i < m_ );
	assert( j >= 0 );	assert( j < n_ );

	top_[ i + j * m_ ] = val;
}

/**
 * Operator overload =
 */
Matrix &Matrix::operator=( const Matrix &T )
{
	assert( m_ == T.m_ );
	assert( n_ == T.n_ );

	for (int i = 0; i < m_ * n_; i++)
		top_[i] = T.top_[i];

	return *this;
}

/**
 * Operator overload []
 *
 * @param i index
 */
double &Matrix::operator[]( const int i ) const
{
	assert( i >= 0 );
	assert( i < m_ * n_ );

	return top_[ i ];
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */
double &Matrix::operator()( const int i, const int j ) const
{
	assert( i >= 0 );  assert( i < m_ );
	assert( j >= 0 );  assert( j < n_ );

	return top_[ i + j * m_ ];
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 */
void Matrix::File_Out( const char* fname )
{
	ofstream matf(fname);
	if (!matf) {
		std::cerr << "Unable to open " << fname << std::endl;
		exit(1);
	}

	matf << m_ << std::endl;
	matf << n_ << std::endl;
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			matf << top_[ i + j * m_ ] << ", ";
		}
		matf << std::endl;
	}
	matf.close();
}



/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 * @param dig number of output digit
 */
void Matrix::File_Out( const char* fname, const unsigned dig )
{
	ofstream matf(fname);
	if (!matf) {
		std::cerr << "Unable to open " << fname << std::endl;
		exit(1);
	}

	matf << m_ << std::endl;
	matf << n_ << std::endl;
	matf.precision(dig);
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			matf << top_[ i + j * m_ ] << " ";
		}
		matf << std::endl;
	}
	matf.close();
}
