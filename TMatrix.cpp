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

	M_ = N_ = mb_ = nb_ = mt_ = nt_ = 0;
	top_ = NULL;
}

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
TMatrix::TMatrix( const int M, const int N,
		const int mb, const int nb,
		const int ib )
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
	mt_ = M_ % mb_ == 0 ? M_ / mb_ : M_ / mb_ + 1;
	nt_ = N_ % nb_ == 0 ? N_ / nb_ : N_ / nb_ + 1;

	try
	{
		top_ = new BMatrix* [ mt_ * nt_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
			int tm = ( i != mt_-1 ) ? mb_ : M_ - i * mb_;
			int tn = ( j != nt_-1 ) ? nb_ : N_ - j * nb_;
			try
			{
				top_[ i + j * mt_ ] = new BMatrix( tm, tn, ib );
			}
			catch (char* eb)
			{
				cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}

			#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
//			cout << "top_[" << i + j * mt_ << "] = " << &top_[ i + j * mt_ ] << endl;
			#endif
		}
}

/**
 * Copy constructor
 *
 * @param T TMatrix object
  */
TMatrix::TMatrix( const TMatrix& T )
{
	#ifdef DEBUG
	cout << "TMatrix(T)\n";
	#endif

	assert( T.top_ != NULL );

	M_ = T.M_;
	N_ = T.N_;
	mb_ = T.mb_;
	nb_ = T.nb_;
	mt_ = M_ % mb_ == 0 ? M_ / mb_ : M_ / mb_ + 1;
	nt_ = N_ % nb_ == 0 ? N_ / nb_ : N_ / nb_ + 1;

	try
	{
		top_ = new BMatrix* [ mt_ * nt_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
		exit(EXIT_FAILURE);
	}

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
			int tm = ( i != mt_-1 ) ? mb_ : M_ - i * mb_;
			int tn = ( j != nt_-1 ) ? nb_ : N_ - j * nb_;
			try
			{
				top_[ i + j * mt_ ] = new BMatrix( tm, tn, ib );
			}
			catch (char* eb)
			{
				cerr << "Can't allocate memory space for TMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}

			#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
//			cout << "top_[" << i + j * mt_ << "] = " << &top_[ i + j * mt_ ] << endl;
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

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
//			cout << "top_[" << i + j * mt_ << "] = " << &top_[ i + j * mt_ ] << endl;
#endif
			delete top_[ i + j * mt_ ];
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
  for (int I = 0; I < N_; I++) {
    for (int J = 0; J < N_; J++) {
      // (ti,tj) : Tile Index
      int ti = I / mb_;
      int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      int i = I % mb_;
      int j = J % nb_;

      top_[ ti + tj * nt_ ]->Set_Val(i, j, Tmp(I,J));
    }
  }
}

/**
 * Operator overload ()
 *
 * @param i low index
 * @param j column index
 */
BMatrix* TMatrix::operator()( const int i, const int j ) const
{
	assert( i > 0 && j > 0 );
	assert( i < mt_ && j < nt_ );

	return top_[ i + j * mt_ ];
}
