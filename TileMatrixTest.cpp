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
	cout << "Matrix block start" << endl;
	{
		Matrix A(100,50);
		cout << "m = " << A.m() << ", n = " << A.n() << endl;
	}
	cout << "Matrix block end" << endl;
	cout << endl << "BMatrix block start" << endl;
	{
		BMatrix A(100,10,5);
		cout << A(0,0) << endl;
		cout << "m = " << A.m() << ", n = " << A.n() << ", ib = " << A.ib() << endl;
	}
	cout << "BMatrix block start" << endl;
	cout << endl << "TMatrix block start" << endl;
	{
		TMatrix A(20,18,10,10,5);
		cout << "M = " << A.M() << ", n = " << A.N() << ", mb = " << A.mb() << ", nb = " << A.nb();
		cout << ", p = " << A.p() << ", q = " << A.q() << endl;

		A.Set_Rnd(20150903);
//		cout << (A[0]) << endl;
	}
	cout << "TMatrix block end" << endl;

	cout << "TileMatrixTest end" << endl;
	return EXIT_SUCCESS;
}
