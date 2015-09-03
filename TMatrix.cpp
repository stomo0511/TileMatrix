/*
 * TMatrix.cpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Matrix.hpp"
#include "TMatrix.hpp"

using namespace std;

/**
 * Default constructor
 */
TMatrix::TMatrix()
{
#ifdef DEBUG
	cout << "TMatrix()\n";
#endif

	M_ = N_ = mb_ = nb_ = p_ = q_ = 0;
	top_ = NULL;
}

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
TMatrix::TMatrix( const unsigned int M, const unsigned int N,
		const unsigned int mb, const unsigned int nb,
		const unsigned int ib )
{
#ifdef DEBUG
	cout << "TMatrix(M,N,mb,nb,ib)\n";
#endif

	assert( M > 0 && N > 0 && mb > 0 && nb > 0 && ib > 0);
	assert( mb < M && nb < N && ib < nb );

	M_ = M;
	N_ = N;
	mb_ = mb;
	nb_ = nb;
	p_ = M_ % mb_ == 0 ? M_ / mb_ : M_ / mb_ + 1;
	q_ = N_ % nb_ == 0 ? N_ / nb_ : N_ / nb_ + 1;

	try
	{
		top_ = new BMatrix* [ p_ * q_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}

	for (unsigned int j=0; j<q_; j++)
		for (unsigned int i=0; i<p_; i++)
		{
			unsigned int tm = ( i != p_-1 ) ? mb_ : M_ - i * mb_;
			unsigned int tn = ( j != q_-1 ) ? nb_ : N_ - j * nb_;
			try
			{
				top_[ i + j * p_ ] = new BMatrix( tm, tn, ib );
			}
			catch (char* eb)
			{
				cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}
#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * p_ ])->m() << ",";
			cout << (top_[ i + j * p_ ])->n() << ",";
			cout << (top_[ i + j * p_ ])->ib() << ")\n";
//			cout << "top_[" << i + j * p_ << "] = " << &top_[ i + j * p_ ] << endl;
#endif
		}
}

/**
 * Destructor
 *
 */
TMatrix::~TMatrix()
{
#ifdef DEBUG
	cout << "\n~TMatrix()\n";
#endif

	for (unsigned int j=0; j<q_; j++)
		for (unsigned int i=0; i<p_; i++)
		{
#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * p_ ])->m() << ",";
			cout << (top_[ i + j * p_ ])->n() << ",";
			cout << (top_[ i + j * p_ ])->ib() << ")\n";
//			cout << "top_[" << i + j * p_ << "] = " << &top_[ i + j * p_ ] << endl;
#endif
			delete top_[ i + j * p_ ];
		}
	delete [] top_;
}

/**
 * Assign random numbers to the elements
 *
 * @param seed Seed of random number generator
 */
void TMatrix::Set_Rnd( const unsigned seed )
{
  Matrix Tmp(M_,N_);
  Tmp.Set_Rnd( seed );

  // (I,J) : Index of the elements of Matrix
  for (unsigned int I = 0; I < N_; I++) {
    for (unsigned int J = 0; J < N_; J++) {
      // (ti,tj) : Tile Index
      unsigned int ti = I / mb_;
      unsigned int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      unsigned int i = I % mb_;
      unsigned int j = J % nb_;

      top_[ ti + tj * q_ ]->Set_Val(i, j, Tmp(I,J));
    }
  }
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */
BMatrix* TMatrix::operator()( const unsigned int i, const unsigned int j ) const
{
	assert( i > 0 && j > 0 );
	assert( i < p_ && j < q_ );

	return top_[ i + j * p_ ];
}
