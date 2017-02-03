/*
 * PivMatrix.cpp
 *
 *  Created on: 2017/02/03
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>

#include "PivMatrix.hpp"

using namespace std;

/**
 * Default constructor
 */
PivMatrix::PivMatrix()
{
	#ifdef DEBUG
	cout << "PivMatrix()\n";
	#endif

	mt_ = nt_ = nb_ = 0;
	top_ = NULL;
}

/**
 * Constructor
 *
 * @param mt number of lows of the pivot matrix
 * @param nt number of columns of the pivot matrix
 * @param nb length of the pivot vector
 */
PivMatrix::PivMatrix( const int mt, const int nt, const int nb )
{
	#ifdef DEBUG
	cout << "PivMatrix(mt,nt,nb)\n";
	#endif

	assert( mt > 0 && nt > 0 && nb > 0 );

	mt_ = mt;
	nt_ = nt;
	nb_ = nb;

	try
	{
		top_ = new int* [ mt_ * nt_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for PivMatrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
			try
			{
				top_[ i + j * mt_ ] = new int [nb_];
			}
			catch (char* eb)
			{
				cerr << "Can't allocate memory space for PivMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}
		}
}

/**
 * Copy constructor
 *
 * @param T TMatrix object
  */
//PivMatrix::PivMatrix( const PivMatrix& T )
//{
//	#ifdef DEBUG
//	cout << "PivMatrix(T)\n";
//	#endif
//
//	assert( T.top_ != NULL );
//
//	mt_ = T.mt_;
//	nt_ = T.nt_;
//	nb_ = T.nb_;
//
//	try
//	{
//		top_ = new int* [ mt_ * nt_ ];
//	}
//	catch (char* eb)
//	{
//		cerr << "Can't allocate memory space for PivMatrix class: " << eb << endl;
//		exit(EXIT_FAILURE);
//	}
//
//	for (int j=0; j<nt_; j++)
//		for (int i=0; i<mt_; i++)
//		{
//			try
//			{
//				top_[ i + j * mt_ ] = new int [nb_];
//			}
//			catch (char* eb)
//			{
//				cerr << "Can't allocate memory space for PivMatrix class: " << eb << endl;
//				exit(EXIT_FAILURE);
//			}
//
//			// Copy all entries of T(i,j)
//			for (int k=0; k<mt_*nt_; k++)
//				(top_[ i + j*mt_ ])->operator [](k) = (T.top_[ i + j*mt_ ])->operator [](k);
//		}
//}

/**
 * Destructor
 *
 */
PivMatrix::~PivMatrix()
{
	#ifdef DEBUG
	cout << "\n~PivMatrix()\n";
	#endif

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
			delete top_[ i + j * mt_ ];
		}
	delete [] top_;
}

/**
 * Operator overload []
 *
 * @param i index
 */
int* PivMatrix::operator[]( const int i ) const
{
	assert( i >= 0 );
	assert( i < mt_ * nt_ );

	return top_[ i ];
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */
int* PivMatrix::operator()( const int i, const int j ) const
{
	assert( i >= 0 && j >= 0 );
	assert( i < mt_ && j < nt_ );

	return top_[ i + j * mt_ ];
}
