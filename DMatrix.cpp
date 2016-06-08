/*
 * TMatrix.cpp
 *
 *  Created on: 2016/06/08
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include "Matrix.hpp"
#include "DMatrix.hpp"

using namespace std;

/**
 * Default constructor
 */
DMatrix::DMatrix()
{
	#ifdef DEBUG
	cout << "DMatrix()\n";
	#endif

	M_ = N_ = P_ = mb_ = nb_ = mt_ = nt_ = 0;
	top_ = NULL;
}

/**
 * Constructor
 *
 * @param M number of lows of the matrix
 * @param N number of columns of the matrix
 * @param P number of domains
 * @param mb number of lows of the tile
 * @param nb number of columns of the tile
 * @param ib width of inner block
 */
DMatrix::DMatrix( const int M, const int N, const int P,
		const int mb, const int nb,
		const int ib )
{
	#ifdef DEBUG
	cout << "DMatrix(M,N,P,mb,nb,ib)\n";
	#endif

	assert( M > 0 && N > 0 && mb > 0 && nb > 0 && ib > 0);
	assert( mb <= M && nb <= N && ib <= nb );
	assert( M % P == 0);
	assert( M / P >= N );

	M_ = M;
	N_ = N;
	P_ = P;
	mb_ = mb;
	nb_ = nb;
	mt_ = M_ % mb_ == 0 ? M_ / mb_ : M_ / mb_ + 1;
	nt_ = N_ % nb_ == 0 ? N_ / nb_ : N_ / nb_ + 1;
	mtl_ = (M_ / P_) / mb_;

	try
	{
		top_ = new BMatrix* [ mt_ * nt_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for DMatrix class: " << eb << endl;
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
				cerr << "Can't allocate memory space for DMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}

			#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
			#endif
		}
}

/**
 * Copy constructor
 *
 * @param T DMatrix object
  */
DMatrix::DMatrix( const DMatrix& T )
{
	#ifdef DEBUG
	cout << "DMatrix(T)\n";
	#endif

	assert( T.top_ != NULL );

	M_ = T.M_;
	N_ = T.N_;
	P_ = T.P_;
	mb_ = T.mb_;
	nb_ = T.nb_;
	mt_ = M_ % mb_ == 0 ? M_ / mb_ : M_ / mb_ + 1;
	nt_ = N_ % nb_ == 0 ? N_ / nb_ : N_ / nb_ + 1;
	mtl_ = M_ / P_;

	int ib = T(0,0,0)->ib();

	try
	{
		top_ = new BMatrix* [ mt_ * nt_ ];
	}
	catch (char* eb)
	{
		cerr << "Can't allocate memory space for DMatrix class: " << eb << endl;
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
				cerr << "Can't allocate memory space for DMatrix class: " << eb << endl;
				exit(EXIT_FAILURE);
			}

			// Copy all entries of T(i,j)
			for (int k=0; k<tm*tn; k++)
				(top_[ i + j*mt_ ])->operator [](k) = (T.top_[ i + j*mt_ ])->operator [](k);

			#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
			#endif
		}
}

//
// Don't work well yet
//
/**
 * Destructor
 *
 */
DMatrix::~DMatrix()
{
	#ifdef DEBUG
	cout << "\n~DMatrix()\n";
	#endif

	for (int j=0; j<nt_; j++)
		for (int i=0; i<mt_; i++)
		{
			#ifdef DEBUG
			cout << "BMatrix(" << (top_[ i + j * mt_ ])->m() << ",";
			cout << (top_[ i + j * mt_ ])->n() << ",";
			cout << (top_[ i + j * mt_ ])->ib() << ")\n";
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
void DMatrix::Set_Rnd( const unsigned seed )
{
	Matrix Tmp(M_,N_);
	Tmp.Set_Rnd( seed );

	// (I,J) : Index of the elements of Matrix
	for (int I = 0; I < M_; I++) {
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
 * Set matrix to the identity matrix
 *
 */
void DMatrix::Set_Iden()
{
	for (int i=0; i<mt_; i++)
		for (int j=0; j<nt_; j++)
			if ( i == j )
				top_[ i + j*mt_ ]->Set_Iden();
			else
				top_[ i + j*mt_ ]->Set_Zero();
}

/**
 * Operator overload ()
 *
 * @param p domain index
 * @param i low index
 * @param j column index
 */
BMatrix* DMatrix::operator()( const int p, const int i, const int j ) const
{
	assert( p >= 0 && p < M_ / P_ );
	assert( i >= 0 && j >= 0 );
	assert( i < mt_ && j < nt_ );

	return top_[ p*mtl_ + i + j * mt_ ];
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 */
void DMatrix::File_Out( const char* fname )
{
	  Matrix Tmp(M_,N_);

	  // (I,J) : Index of the elements of Matrix
	  for (int I = 0; I < M_; I++) {
		for (int J = 0; J < N_; J++) {
		  // (ti,tj) : Tile Index
		  int ti = I / mb_;
		  int tj = J / nb_;
		  // (i,j) : Index of the elements of Tile
		  int i = I % mb_;
		  int j = J % nb_;

		  double val = top_[ ti + tj*mt_ ]->operator()(i, j);
		  Tmp.Set_Val( I, J, val );
		}
	  }
	  Tmp.File_Out( fname );
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 * @param dig number of output digit
 */
void DMatrix::File_Out( const char* fname, const unsigned dig )
{
	Matrix Tmp(M_,N_);

	// (I,J) : Index of the elements of Matrix
	for (int I = 0; I < M_; I++) {
		for (int J = 0; J < N_; J++) {
			// (ti,tj) : Tile Index
			int ti = I / mb_;
			int tj = J / nb_;
			// (i,j) : Index of the elements of Tile
			int i = I % mb_;
			int j = J % nb_;

			double val = top_[ ti + tj*mt_ ]->operator()(i, j);
			Tmp.Set_Val( I, J, val );
		}
	}
	Tmp.File_Out( fname, dig );
}

/*
 * copy matrix elements to the Matrix class
 */
void DMatrix::Mat_Copy( Matrix& A )
{
	assert( M_ == A.m() ); assert( N_ == A.n() );

	double val;
	for (int I=0; I<M_; I++ )
		for (int J=0; J<N_; J++)
		{
			// (ti,tj) : Tile index
			int ti = I / mb_;
			int tj = J / nb_;

			// (i,j) : Index of the tile elements
			int i = I % mb_;
			int j = J % nb_;

			val = top_[ ti + tj*mt_ ]->operator ()(i,j);
			A.Set_Val(I,J,val);
		}
}

/*
 * copy matrix elements to the array
 */
void DMatrix::Array_Copy( double* array )
{
	assert( array != NULL );

	for (int I=0; I<M_; I++ )
		for (int J=0; J<N_; J++)
		{
			// (ti,tj) : Tile index
			int ti = I / mb_;
			int tj = J / nb_;

			// (i,j) : Index of the tile elements
			int i = I % mb_;
			int j = J % nb_;

			array[ I + J*M_ ] = top_[ ti + tj*mt_ ]->operator ()(i,j);
		}
}

