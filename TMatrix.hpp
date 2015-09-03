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
	unsigned int M_;	// number of lows of the matrix
	unsigned int N_;	// number of columns of the matrix
	unsigned int mb_;	// number of lows of the tile
	unsigned int nb_;	// number of columns of the tile
	unsigned int p_;    // number of low tiles
	unsigned int q_;    // number of column tiles

public:
	// Default constructor
	TMatrix();

	// Constructor
	TMatrix( const unsigned int M, const unsigned int N,
			const unsigned int mb, const unsigned int nb,
			const unsigned int ib );

	// Destructor
	virtual ~TMatrix();

	// Getters
	unsigned int M()  const { return M_; }
	unsigned int N()  const { return N_; }
	unsigned int mb() const { return mb_; }
	unsigned int nb() const { return nb_; }
	unsigned int p()  const { return p_; }
	unsigned int q()  const { return q_; }

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Operator overload
	BMatrix* operator[]( const unsigned int i ) const;
	BMatrix* operator()( const unsigned int i, const unsigned int j ) const;
};

#endif /* TMATRIX_HPP_ */
