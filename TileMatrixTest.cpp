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
#include "DMatrix.hpp"

using namespace std;

int main( int argc, char* argv[] )
{
//	cout << "Matrix block start" << endl;
//	{
//		Matrix A(25,25);
//		cout << "m = " << A.m() << ", n = " << A.n() << endl;
//		A.Set_Rnd(20151121);
//		A.Show_all();
//		A.Set_Iden();
//		A.Show_all();
//		A.Set_Zero();
//		A.Show_all();
//
//		Matrix B = A;
//		B.Show_all();
//	}
//	cout << "Matrix block end" << endl;

/*	cout << endl << "<<< BMatrix block start >>>" << endl;
	{
		BMatrix A(10,10,5);
		cout << A(0,0) << endl;
		cout << "m = " << A.m() << ", n = " << A.n() << ", ib = " << A.ib() << endl;
		A.Set_Rnd(20151121);
		A.Show_all();

		A.Set_Iden();
		A.Show_all();

		A.Set_Zero();
		A.Show_all();

		A.Set_Rnd(20151121);
		BMatrix B = A;
		B.Show_all();
		cout << "B.m = " << B.m() << ", B.n = " << B.n() << ", B.ib = " << B.ib() << endl;
	}
	cout << "<<< BMatrix block end >>>" << endl;*/

	cout << endl << "<<< TMatrix block start >>>" << endl;
	{
		TMatrix A(20,10,10,10,5);
		cout << "M = " << A.M() << ", n = " << A.N() << ", mb = " << A.mb() << ", nb = " << A.nb();
		cout << ", mt = " << A.mt() << ", nt = " << A.nt() << endl << endl;

		A.Set_Rnd(20151121);
		A(0,0)->Show_all();
		A(1,0)->Show_all();

//		A.Set_Iden();
//		A(1,1)->Show_all();

		A.Set_Rnd(20151121);
		TMatrix B = A;
		B(0,0)->Show_all();
	}
	cout << "<<< TMatrix block end >>>" << endl;

	cout << endl << "<<< DMatrix block start >>>" << endl;
	{
		DMatrix A(20,10,2,10,10,5);
		cout << "M = " << A.M() << ", n = " << A.N() << ", mb = " << A.mb() << ", nb = " << A.nb();
		cout << ", P = " << A.P() << ", mtl = " << A.mtl();
		cout << ", mt = " << A.mt() << ", nt = " << A.nt() << endl << endl;

		A.Set_Rnd(20151121);
		A(0,0,0)->Show_all();
		A(1,0,0)->Show_all();

		A.Set_Iden();
		A(0,0,0)->Show_all();

		A.Set_Rnd(20151121);
		DMatrix B = A;
		B(0,0,0)->Show_all();
	}
	cout << "<<< DMatrix block end >>>" << endl;

	cout << "TileMatrixTest end" << endl;
	return EXIT_SUCCESS;
}
