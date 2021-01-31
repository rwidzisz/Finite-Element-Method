#ifndef HEAT_H
#define HEAT_H

#include "Heat.h"

using namespace std;

double *heat::jakobian_C(double n, double E, element_siatki element) {

	double N1 = 0.25*(1 - E)*(1 - n);
	double N2 = 0.25*(1 + E)*(1 - n);
	double N3 = 0.25*(1 + E)*(1 + n);
	double N4 = 0.25*(1 - E)*(1 + n);

	double dN1dE = -0.25*(1 - n);
	double dN2dE = 0.25*(1 - n);
	double dN3dE = 0.25*(1 + n);
	double dN4dE = -0.25*(1 + n);

	double dN1dn = -0.25*(1 - E);
	double dN2dn = -0.25*(1 + E);
	double dN3dn = 0.25*(1 + E);
	double dN4dn = 0.25*(1 - E);

	double dxdE = dN1dE * element.pkt1.x + dN2dE * element.pkt2.x + dN3dE * element.pkt3.x + dN4dE * element.pkt4.x;
	double dydn = dN1dn * element.pkt1.y + dN2dn * element.pkt2.y + dN3dn * element.pkt3.y + dN4dn * element.pkt4.y;
	double dydE = dN1dE * element.pkt1.y + dN2dE * element.pkt2.y + dN3dE * element.pkt3.y + dN4dE * element.pkt4.y;
	double dxdn = dN1dn * element.pkt1.x + dN2dn * element.pkt2.x + dN3dn * element.pkt3.x + dN4dn * element.pkt4.x;


	//Obliczanie jakobianu
	double macierz_jakobiego[4] = { dxdE ,dxdn ,dydE ,dydn };
	double det = (macierz_jakobiego[0] * macierz_jakobiego[3]) - ((macierz_jakobiego[1] * macierz_jakobiego[2]));

	//cout <<endl<< N1 << "\t" << N2 << "\t" << N3 << "\t" << N4 << endl;


	double *tab = new double[5];
	tab[0] = N1;
	tab[1] = N2;
	tab[2] = N3;
	tab[3] = N4;
	tab[4] = det;

	return tab;
}

double **heat::macierzC(double tab[5], int cp) //tablica z danymi z liczenia jakobianu
{
	double N[1][4];
	double NT[4][1];
	double det = tab[4];

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

	//Tworzenie macierzy C
	double** matrixC = 0;
	matrixC = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixC[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixC[i][j] = matrixN[i][j] * det*cp;
		}

	}

	return matrixC;
}

double **heat::macierz_lokalna_C_dwupunkt(vector<pcalkowania> ZPC, element_siatki element, int cp, int ro) {

	vector<pcalkowania>zbior_pcalkowania = ZPC;

	double matrixCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixCL[i][j] = 0;
		}
	}

	for (int pc_integer = 0; pc_integer < 4; pc_integer++)
	{
		double *tab_dynamic = jakobian_C(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabC = macierzC(tab_dynamic,cp);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixCL[i][j] = matrixCL[i][j] + (tabC[i][j] * ro);
			}
		}

	}

/*	cout << "\n\nLokalna Macierz C" <<endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixCret = 0;
	matrixCret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixCret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixCret[i][j] = matrixCL[i][j];
		}

	}


	return matrixCret;
}

double **heat::macierz_lokalna_C_trzypunkt(vector<pcalkowania> ZPC, element_siatki element, int cp, int ro) {

	vector<pcalkowania>zbior_pcalkowania = ZPC;

	vector<double>zbior_wag;

	double wg[3] = { (5.0 / 9.0),	(8.0 / 9.0), (5.0 / 9.0) };
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double current_waga = wg[i] * wg[j];
			zbior_wag.push_back(current_waga);
		}
	}

	double matrixCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixCL[i][j] = 0;

		}
	}

	for (int pc_integer = 0; pc_integer < 9; pc_integer++) //9 punkotw ca³kowania
	{
		double *tab_dynamic = jakobian_C(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabC = macierzC(tab_dynamic,cp);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixCL[i][j] = matrixCL[i][j] + (tabC[i][j] * ro * zbior_wag[pc_integer]);
			}
		}

	}

	/*cout << "\n\nLokalna Macierz C" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixCret = 0;
	matrixCret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixCret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixCret[i][j] = matrixCL[i][j];
		}

	}

	return matrixCret;
}

double **heat::macierz_lokalna_C_czteropunkt(vector<pcalkowania> ZPC, element_siatki element, int cp,int ro) {

	vector<pcalkowania>zbior_pcalkowania = ZPC;

	vector<double>zbior_wag;
	double wag_value1 = (18 + sqrt(30)) / 36;
	double wag_value2 = (18 - sqrt(30)) / 36;
	double wg[4] = { wag_value2,wag_value1,wag_value1,wag_value2 };
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 4; j++) {
			double current_waga = wg[i] * wg[j];
			zbior_wag.push_back(current_waga);
		}
	}

	double matrixCL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixCL[i][j] = 0;

		}
	}

	for (int pc_integer = 0; pc_integer < 16; pc_integer++) //16 punkotw ca³kowania
	{
		double *tab_dynamic = jakobian_C(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabC = macierzC(tab_dynamic,cp);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixCL[i][j] = matrixCL[i][j] + (tabC[i][j] * ro * zbior_wag[pc_integer]);
			}
		}

	}

	/*cout << "\n\nLokalna Macierz C" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixCL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixCret = 0;
	matrixCret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixCret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixCret[i][j] = matrixCL[i][j];
		}

	}

	return matrixCret;
}


#endif