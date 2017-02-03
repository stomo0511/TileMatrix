/*
 * PivMatrix.hpp
 *
 *  Created on: 2017/02/03
 *      Author: stomo
 */

#ifndef PIVMATRIX_HPP_
#define PIVMATRIX_HPP_

class PivMatrix {

protected:
	int** top_;		// pointer to the pivot matrix
	int mt_;		// number of lows
	int nt_;		// number of columns
	int nb_;		// pivot vector length

public:
	// Default constructor
	PivMatrix();

	// Constructor
	PivMatrix( const int mt, const int nt, const int nb );

	// Copy constructor
	//PivMatrix( const PivMatrix& T );

	// Destructor
	virtual ~PivMatrix();

	// Getters
	int mt() const { return mt_; }
	int nt() const { return nt_; }
	int nb() const { return nb_; }

	// Operator overload
	int* operator[]( const int i ) const;
	int* operator()( const int i, const int j ) const;
};


#endif /* PIVMATRIX_HPP_ */
