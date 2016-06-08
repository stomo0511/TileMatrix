/*
 * DMatrix.hpp
 *
 *  Created on: 2016/06/08
 *      Author: stomo
 */

#ifndef DMATRIX_HPP_
#define DMATRIX_HPP_

#include "BMatrix.hpp"

class DMatrix {

private:
	BMatrix** top_;	// pointer to the BMatrix
	int M_;			// number of lows of the matrix
	int N_;			// number of columns of the matrix
	int P_;			// number of domains
	int mb_;		// number of lows of the tile
	int nb_;		// number of columns of the tile
	int mt_;    	// number of low tiles
	int nt_;    	// number of column tiles
	int mtl_;		// number of low tiles of the domain

public:
	// Default constructor
	DMatrix();

	// Constructor
	DMatrix( const int M, const int N, const int P,
			const int mb, const int nb,
			const int ib );

	// Copy constructor
	DMatrix( const DMatrix& T );

	// Destructor
	virtual ~DMatrix();

	// Getters
	int M()  const { return M_; }
	int N()  const { return N_; }
	int P()  const { return P_; }
	int mb() const { return mb_; }
	int nb() const { return nb_; }
	int mt() const { return mt_; }
	int nt() const { return nt_; }
	int mtl() const { return mtl_; }

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Set matrix to the identity matrix
	void Set_Iden();

	// Operator overload
	BMatrix* operator[]( const int i ) const;
	BMatrix* operator()( const int p, const int i, const int j ) const;

	// Save DMatrix to the file
	void File_Out( const char* fname );  // Not Yet
	void File_Out( const char* fname, const unsigned dig );  // Not Yet

	// copy matrix elements to the Matrix class
	void Mat_Copy( Matrix& A );

	// copy matrix elements to the array
	void Array_Copy( double *array );
};

#endif /* DMATRIX_HPP_ */
