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
	BMatrix** top_;		// pointer to the BMatrix
	int M_;	// number of lows of the matrix
	int N_;	// number of columns of the matrix
	int mb_;	// number of lows of the tile
	int nb_;	// number of columns of the tile
	int p_;    // number of low tiles
	int q_;    // number of column tiles

public:
	// Default constructor
	TMatrix();

	// Constructor
	TMatrix( const int M, const int N,
			const int mb, const int nb,
			const int ib );

	// Destructor
	virtual ~TMatrix();

	// Getters
	int M()  const { return M_; }
	int N()  const { return N_; }
	int mb() const { return mb_; }
	int nb() const { return nb_; }
	int p()  const { return p_; }
	int q()  const { return q_; }

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Operator overload
	BMatrix* operator[]( const int i ) const;
	BMatrix* operator()( const int i, const int j ) const;
};

#endif /* TMATRIX_HPP_ */
