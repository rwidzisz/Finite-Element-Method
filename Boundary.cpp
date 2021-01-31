#ifndef BOUNDARY_H
#define BOUNDARY_H

#include "Boundary.h"

using namespace std;

double **boundary_condition::macierzHBC(double N1, double N2, double N3, double N4) {

	double tab[4] = { N1,N2,N3,N4 };
	double N[1][4];
	double NT[4][1];

	for (int i = 0; i < 4; i++)
	{
		N[0][i] = tab[i];		//rows1 columns4
		NT[i][0] = tab[i];		//rows4 columns1
	}

	double matrixN[4][4];

	//Tworzenie macierzy 
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			matrixN[i][j] = 0;
		}

	//Mno¿enie wektorów, wype³nianie macierzy
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			for (int n = 0; n < 1; ++n)
			{
				matrixN[i][j] += NT[i][n] * N[n][j];
			}
		}
	}

	//Tworzenie macierzy HBC
	double** matrixHBC = 0;
	matrixHBC = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixHBC[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHBC[i][j] = matrixN[i][j];
		}

	}

	return matrixHBC;
}

double **boundary_condition::macierz_lokalna_BC_dwupunkt(element_siatki element, int alfa) {

	double p057 = sqrt(1.0 / 3.0);	//0.57
	double n057 = -1 * p057;
	pcalkowania p1(n057, -1), p2(p057, -1), p3(1, n057), p4(1, p057), p5(0.57, 1), p6(n057, 1), p7(-1, p057), p8(-1, n057);


	double matrixHBCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHBCL[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N11 = 0.5*(1 - p1.x);
		double N21 = 0.5*(1 + p1.x);
		double N31 = 0;
		double N41 = 0;

		double N12 = 0.5*(1 - p2.x);
		double N22 = 0.5*(1 + p2.x);
		double N32 = 0;
		double N42 = 0;

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] + tabH2[i][j])*det*alfa;
			}
		}

	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N11 = 0;
		double N21 = 0.5*(1 - p3.y);
		double N31 = 0.5*(1 + p3.y);
		double N41 = 0;

		double N12 = 0;
		double N22 = 0.5*(1 - p4.y);
		double N32 = 0.5*(1 + p4.y);
		double N42 = 0;

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] + tabH2[i][j])*det* alfa;
			}
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{

		double N11 = 0;
		double N21 = 0;
		double N31 = 0.5*(1 - p5.x);
		double N41 = 0.5*(1 + p5.x);

		double N12 = 0;
		double N22 = 0;
		double N32 = 0.5*(1 - p6.x);
		double N42 = 0.5*(1 + p6.x);

		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] + tabH2[i][j])*det* alfa;
			}
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true)
	{

		double N11 = 0.5*(1 + p7.y);
		double N21 = 0;
		double N31 = 0;
		double N41 = 0.5*(1 - p7.y);

		double N12 = 0.5*(1 + p8.y);
		double N22 = 0;
		double N32 = 0;
		double N42 = 0.5*(1 - p8.y);


		double det = (element.pkt4.y - element.pkt1.y) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] + tabH2[i][j])*det* alfa;
			}
		}
	}

	/*cout << "\n\nLokalna Macierz HBC"<< endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHBCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/



	double** matrixHret = 0;
	matrixHret = new double*[4];
	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHBCL[i][j];
		}

	}
	

	return matrixHret;
}

double **boundary_condition::macierz_lokalna_BC_trzypunkt(element_siatki element, int alfa) {

	double p077 = sqrt(3.0 / 5.0);	//0.77
	double n077 = -1 * p077;	//-0.77
	pcalkowania p1(n077, -1), p2(0, -1), p3(p077, -1),
				p4(-1, n077), p5(-1, 0), p6(-1, p077),
				p7(p077, 1), p8(0, 1), p9(n077, 1),
				p10(-1, p077), p11(-1, 0), p12(-1, n077);


	double wg[3] = { (5.0 / 9.0),	(8.0 / 9.0), (5.0 / 9.0) };


	double matrixHBCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHBCL[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N11 = 0.5*(1 - p1.x);
		double N21 = 0.5*(1 + p1.x);
		double N31 = 0;
		double N41 = 0;

		double N12 = 0.5*(1 - p2.x);
		double N22 = 0.5*(1 + p2.x);
		double N32 = 0;
		double N42 = 0;

		double N13 = 0.5*(1 - p3.x);
		double N23 = 0.5*(1 + p3.x);
		double N33 = 0;
		double N43 = 0;

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2])*det*alfa;
			}
		}

	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N11 = 0;
		double N21 = 0.5*(1 - p4.y);
		double N31 = 0.5*(1 + p4.y);
		double N41 = 0;

		double N12 = 0;
		double N22 = 0.5*(1 - p5.y);
		double N32 = 0.5*(1 + p5.y);
		double N42 = 0;

		double N13 = 0;
		double N23 = 0.5*(1 - p6.y);
		double N33 = 0.5*(1 + p6.y);
		double N43 = 0;

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;


		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2])*det*alfa;
			}
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{

		double N11 = 0;
		double N21 = 0;
		double N31 = 0.5*(1 - p7.x);
		double N41 = 0.5*(1 + p7.x);

		double N12 = 0;
		double N22 = 0;
		double N32 = 0.5*(1 - p8.x);
		double N42 = 0.5*(1 + p8.x);

		double N13 = 0;
		double N23 = 0;
		double N33 = 0.5*(1 - p9.x);
		double N43 = 0.5*(1 + p9.x);

		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2])*det*alfa;
			}
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true)
	{

		double N11 = 0.5*(1 + p10.y);
		double N21 = 0;
		double N31 = 0;
		double N41 = 0.5*(1 - p10.y);

		double N12 = 0.5*(1 + p11.y);
		double N22 = 0;
		double N32 = 0;
		double N42 = 0.5*(1 - p11.y);

		double N13 = 0.5*(1 + p12.y);
		double N23 = 0;
		double N33 = 0;
		double N43 = 0.5*(1 - p12.y);

		double det = (element.pkt4.y - element.pkt1.y) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2])*det*alfa;
			}
		}
	}

	/*cout << "\n\nLokalna Macierz HBC" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHBCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/



	double** matrixHret = 0;
	matrixHret = new double*[4];
	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHBCL[i][j];
		}

	}


	return matrixHret;
}

double **boundary_condition::macierz_lokalna_BC_czteropunkt(element_siatki element, int alfa) {

	double p086 = sqrt((3.0 / 7.0 + (2.0 / 7.0*sqrt(6.0 / 5.0))));
	double n086 = -1 * p086;
	double p033 = sqrt((3.0 / 7.0 - (2.0 / 7.0*sqrt(6.0 / 5.0))));
	double n033 = -1 * p033;
	pcalkowania p1(n086, -1), p2(n033, -1), p3(p033, -1), p4(p086, -1),
				p5(-1, n086), p6(-1, n086), p7(-1, p033), p8(-1, p086),
				p9(p086, 1), p10(p033, 1), p11(n033, 1), p12(n086, 1),
				p13(1, p086), p14(1, p033), p15(1, n033), p16(1, n086);


	double wag_value1 = (18 + sqrt(30)) / 36;
	double wag_value2 = (18 - sqrt(30)) / 36;
	double wg[4] = { wag_value2,wag_value1,wag_value1,wag_value2 };


	double matrixHBCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHBCL[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N11 = 0.5*(1 - p1.x);
		double N21 = 0.5*(1 + p1.x);
		double N31 = 0;
		double N41 = 0;

		double N12 = 0.5*(1 - p2.x);
		double N22 = 0.5*(1 + p2.x);
		double N32 = 0;
		double N42 = 0;

		double N13 = 0.5*(1 - p3.x);
		double N23 = 0.5*(1 + p3.x);
		double N33 = 0;
		double N43 = 0;

		double N14 = 0.5*(1 - p4.x);
		double N24 = 0.5*(1 + p4.x);
		double N34 = 0;
		double N44 = 0;

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);
		double **tabH4 = macierzHBC(N14, N24, N34, N44);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2] + tabH4[i][j] * wg[3])*det*alfa;
			}
		}

	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N11 = 0;
		double N21 = 0.5*(1 - p5.y);
		double N31 = 0.5*(1 + p5.y);
		double N41 = 0;

		double N12 = 0;
		double N22 = 0.5*(1 - p6.y);
		double N32 = 0.5*(1 + p6.y);
		double N42 = 0;

		double N13 = 0;
		double N23 = 0.5*(1 - p7.y);
		double N33 = 0.5*(1 + p7.y);
		double N43 = 0;

		double N14 = 0;
		double N24 = 0.5*(1 - p8.y);
		double N34 = 0.5*(1 + p8.y);
		double N44 = 0;

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;


		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);
		double **tabH4 = macierzHBC(N14, N24, N34, N44);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2] + tabH4[i][j] * wg[3])*det*alfa;
			}
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{

		double N11 = 0;
		double N21 = 0;
		double N31 = 0.5*(1 - p9.x);
		double N41 = 0.5*(1 + p9.x);

		double N12 = 0;
		double N22 = 0;
		double N32 = 0.5*(1 - p10.x);
		double N42 = 0.5*(1 + p10.x);

		double N13 = 0;
		double N23 = 0;
		double N33 = 0.5*(1 - p11.x);
		double N43 = 0.5*(1 + p11.x);

		double N14 = 0;
		double N24 = 0;
		double N34 = 0.5*(1 - p12.x);
		double N44 = 0.5*(1 + p12.x);


		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);
		double **tabH4 = macierzHBC(N14, N24, N34, N44);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2] + tabH4[i][j] * wg[3])*det*alfa;
			}
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true) {

		double N11 = 0.5*(1 + p13.y);
		double N21 = 0;
		double N31 = 0;
		double N41 = 0.5*(1 - p13.y);

		double N12 = 0.5*(1 + p14.y);
		double N22 = 0;
		double N32 = 0;
		double N42 = 0.5*(1 - p14.y);

		double N13 = 0.5*(1 + p15.y);
		double N23 = 0;
		double N33 = 0;
		double N43 = 0.5*(1 - p15.y);

		double N14 = 0.5*(1 + p16.y);
		double N24 = 0;
		double N34 = 0;
		double N44 = 0.5*(1 - p16.y);

		double det = (element.pkt4.y - element.pkt1.y) / 2.0;

		double **tabH1 = macierzHBC(N11, N21, N31, N41);
		double **tabH2 = macierzHBC(N12, N22, N32, N42);
		double **tabH3 = macierzHBC(N13, N23, N33, N43);
		double **tabH4 = macierzHBC(N14, N24, N34, N44);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHBCL[i][j] = matrixHBCL[i][j] + (tabH1[i][j] * wg[0] + tabH2[i][j] * wg[1] + tabH3[i][j] * wg[2] + tabH4[i][j] * wg[3])*det*alfa;
			}
		}
	}

	/*cout << "\n\nLokalna Macierz HBC" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHBCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/



	double** matrixHret = 0;
	matrixHret = new double*[4];
	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHBCL[i][j];
		}

	}

	return matrixHret;
}


#endif