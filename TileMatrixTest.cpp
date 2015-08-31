/*
 * TileMatrixTest.cpp
 *
 *  Created on: 2015/08/31
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include "Matrix.hpp"
#include "BMatrix.hpp"
#include "TMatrix.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
	{
		Matrix A(10,10);
		cout << A(0,0) << endl;
		cout << "Matrix block end" << endl;
	}
	{
		Matrix A(100,50);
		cout << "m = " << A.m() << ", n = " << A.n() << endl;
		cout << "Matrix block end" << endl;
	}
	{
		BMatrix A(100,10,5);
		cout << A(0,0) << endl;
		cout << "m = " << A.m() << ", n = " << A.n() << ", ib = " << A.ib() << endl;
		cout << "BMatrix block end" << endl;
	}
	for (unsigned int i=0; i<5; i++)
	{
		TMatrix A(30000,30000,500,500,5);
		cout << "M = " << A.M() << ", n = " << A.N() << ", mb = " << A.mb() << ", nb = " << A.nb();
		cout << ", p = " << A.p() << ", q = " << A.q() << endl;
		cout << "TMatrix block end" << endl;
	}

	cout << "TileMatrixTest end" << endl;
	return EXIT_SUCCESS;
}
