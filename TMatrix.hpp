/*
 * TMatrix.hpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#ifndef TMATRIX_HPP_
#define TMATRIX_HPP_

#include "BMatrix.hpp"

class TMatrix {

private:
	BMatrix** top_;	// pointer to the BMatrix
	int M_;			// number of lows of the matrix
	int N_;			// number of columns of the matrix
	int mb_;		// number of lows of the tile
	int nb_;		// number of columns of the tile
	int mt_;    	// number of low tiles
	int nt_;    	// number of column tiles

public:
	// Default constructor
	TMatrix();

	// Constructor
	TMatrix( const int M, const int N,
			const int mb, const int nb,
			const int ib );

	// Copy constructor
	TMatrix( const TMatrix& T );

	// Destructor
	virtual ~TMatrix();

	// Getters
	int M()  const { return M_; }
	int N()  const { return N_; }
	int mb() const { return mb_; }
	int nb() const { return nb_; }
	int mt() const { return mt_; }
	int nt() const { return nt_; }

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Set matrix to the identity matrix
	void Set_Iden();

	// Operator overload
	BMatrix* operator[]( const int i ) const;
	BMatrix* operator()( const int i, const int j ) const;

	// Save TMatrix to the file
	void File_Out( const char* fname );  // Not Yet
	void File_Out( const char* fname, const unsigned dig );  // Not Yet

	// copy matrix elements to the Matrix class
	void Mat_Copy( Matrix& A );

	// copy matrix elements to the array
	void Array_Copy( double *array );
};

#endif /* TMATRIX_HPP_ */
