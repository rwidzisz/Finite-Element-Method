#ifndef INTEGRAL_H
#define INTEGRAL_H

#include "Integral.h"

using namespace std;
using std::vector;

double *integral::jakobian(double n, double E, element_siatki element) {


	double dN1dE = -0.25*(1 - n);
	double dN2dE = 0.25*(1 - n);
	double dN3dE = 0.25*(1 + n);
	double dN4dE = -0.25*(1 + n);

	double dN1dn = -0.25*(1 - E);
	double dN2dn = -0.25*(1 + E);
	double dN3dn = 0.25*(1 + E);
	double dN4dn = 0.25*(1 - E);

	//cout << n << " " << E << ": " << dN1dE << " " << dN2dE << " " << dN3dE << " " << dN4dE << endl;

	double dxdE = dN1dE * element.pkt1.x + dN2dE * element.pkt2.x + dN3dE * element.pkt3.x + dN4dE * element.pkt4.x;
	double dydn = dN1dn * element.pkt1.y + dN2dn * element.pkt2.y + dN3dn * element.pkt3.y + dN4dn * element.pkt4.y;
	double dydE = dN1dE * element.pkt1.y + dN2dE * element.pkt2.y + dN3dE * element.pkt3.y + dN4dE * element.pkt4.y;
	double dxdn = dN1dn * element.pkt1.x + dN2dn * element.pkt2.x + dN3dn * element.pkt3.x + dN4dn * element.pkt4.x;


	//Obliczanie jakobianu
	double macierz_jakobiego[4] = { dxdE ,dxdn ,dydE ,dydn };
	double det = (macierz_jakobiego[0] * macierz_jakobiego[3]) - ((macierz_jakobiego[1] * macierz_jakobiego[2]));
	double odwrocona_macierz_jakobiego[4] = { dydn,dydE,dxdn ,dxdE };

	
	double dN1dx = (1 / det) * (dN1dE*dydn + dN1dn * (-dydE));
	double dN2dx = (1 / det) * (dN2dE*dydn + dN2dn * (-dydE));
	double dN3dx = (1 / det) * (dN3dE*dydn + dN3dn * (-dydE));
	double dN4dx = (1 / det) * (dN4dE*dydn + dN4dn * (-dydE));

	double dN1dy = (1 / det) * (dN1dE*(-dxdn) + dN1dn * (dxdE));
	double dN2dy = (1 / det) * (dN2dE*(-dxdn) + dN2dn * (dxdE));
	double dN3dy = (1 / det) * (dN3dE*(-dxdn) + dN3dn * (dxdE));
	double dN4dy = (1 / det) * (dN4dE*(-dxdn) + dN4dn * (dxdE));

	//Drukowanie danych

/*
	cout << endl << "\n\nMacierz jakobiego:" << endl;
	for (int i = 0; i < 4; i++) {
		cout << macierz_jakobiego[i] << "\t";
		if (i == 1) {
			cout << endl;
		}
	}
	cout << endl << "\nWyznacznik:" << det << endl;
	cout << "\nOdwrocona macierz jakobiego:" << endl;
	for (int i = 0; i < 4; i++) {
		cout << odwrocona_macierz_jakobiego[i] << "\t";
		if (i == 1) {
			cout << endl;
		}
	}

	cout << "\n \n Wyniki dla punktu calkowania" << endl;
	cout << " dN1dx: " << dN1dx << "\t" << " dN2dx: " << dN2dx << "\t" << " dN3dx: " << dN3dx << "\t" << " dN4dx: " << dN4dx << endl;
	cout << " dN1dy: " << dN2dy << "\t" << " dN2dy: " << dN2dy << "\t" << " dN3dy: " << dN3dy << "\t" << " dN4dy: " << dN4dy << endl;
*/


	double *tab = new double[9];
	tab[0] = dN1dx;
	tab[1] = dN2dx;
	tab[2] = dN3dx;
	tab[3] = dN4dx;
	tab[4] = dN1dy;
	tab[5] = dN2dy;
	tab[6] = dN3dy;
	tab[7] = dN4dy;
	tab[8] = det;

	return tab;
}

double **integral::macierzH(double tab[9],int k) //tablica z danymi z liczenia jakobianu
{
	double dNPx[1][4];
	double dNPxT[4][1];
	double dNPy[1][4];
	double dNPyT[4][1];
	double det = tab[8];

	for (int i = 0; i < 4; i++)
	{
		dNPx[0][i] = tab[i];		//rows1 columns4
		dNPxT[i][0] = tab[i];		//rows4 columns1

		dNPy[0][i] = tab[i + 4];	//rows1 columns4
		dNPyT[i][0] = tab[i + 4];	//rows4 columns1
	}

	double matrixX[4][4], matrixY[4][4];

	//Tworzenie macierzy 
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			matrixX[i][j] = 0;
			matrixY[i][j] = 0;
		}

	//Mno¿enie wektorów, wype³nianie macierzy
	for (int i = 0; i < 4; ++i)
	{
		for (int j = 0; j < 4; ++j)
		{
			for (int n = 0; n < 1; ++n)
			{
				matrixX[i][j] += dNPxT[i][n] * dNPx[n][j];
				matrixY[i][j] += dNPyT[i][n] * dNPy[n][j];
			}
		}
	}

	//Tworzenie macierzy H
	double** matrixH = 0;
	matrixH = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixH[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixH[i][j] = (matrixX[i][j] + matrixY[i][j])*det*k;
		}

	}

	return matrixH;
}

double **integral::macierz_lokalna_dwupunkt(vector<pcalkowania> ZPC, element_siatki element, int k) {

	vector<pcalkowania>zbior_pcalkowania = ZPC;

	double matrixHL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHL[i][j] = 0;
		}
	}

	for (int pc_integer = 0; pc_integer < 4; pc_integer++)
	{
		double *tab_dynamic = jakobian(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabH = macierzH(tab_dynamic,k);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHL[i][j] = matrixHL[i][j] + tabH[i][j];
			}
		}

	}


	/*cout << "\n\nLokalna Macierz H:" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixHret = 0;
	matrixHret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHL[i][j];
		}

	}

	
	return matrixHret;
}

double **integral::macierz_lokalna_trzypunkt(vector<pcalkowania> ZPC, element_siatki element, int k) {

	vector<pcalkowania>zbior_pcalkowania = ZPC;

	vector<double>zbior_wag;

	double wg[3] = { (5.0 / 9.0),	(8.0 / 9.0), (5.0 / 9.0) };
	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++) {
			double current_waga = wg[i] * wg[j];
			zbior_wag.push_back(current_waga);
		}
	}

	double matrixHL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHL[i][j] = 0;

		}
	}

	for (int pc_integer = 0; pc_integer < 9; pc_integer++) //9 punkotw ca³kowania
	{
		double *tab_dynamic = jakobian(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabH = macierzH(tab_dynamic,k);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHL[i][j] = matrixHL[i][j] + (tabH[i][j] * zbior_wag[pc_integer]);
			}
		}

	}

	/*cout << "\n\nLokalna Macierz H:" << endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixHret = 0;
	matrixHret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHL[i][j];
		}

	}

	return matrixHret;
}

double **integral::macierz_lokalna_czteropunkt(vector<pcalkowania> ZPC, element_siatki element, int k) {

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

	double matrixHL[4][4];
	for (int i = 0; i < 4; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixHL[i][j] = 0;

		}
	}

	for (int pc_integer = 0; pc_integer < 16; pc_integer++) //16 punkotw ca³kowania
	{
		double *tab_dynamic = jakobian(zbior_pcalkowania[pc_integer].x, zbior_pcalkowania[pc_integer].y, element);
		double **tabH = macierzH(tab_dynamic,k);

		for (int i = 0; i < 4; ++i) {
			for (int j = 0; j < 4; ++j)
			{
				matrixHL[i][j] = matrixHL[i][j] + (tabH[i][j] * zbior_wag[pc_integer]);
			}
		}

	}

	/*cout << "\n\nLokalna Macierz H:"<< endl;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixHL[i][j];
			if (j == 4 - 1)
				cout << endl;
		}*/

	double** matrixHret = 0;
	matrixHret = new double*[4];

	for (int i = 0; i < 4; ++i) {
		matrixHret[i] = new double[4];

		for (int j = 0; j < 4; ++j)
		{
			matrixHret[i][j] = matrixHL[i][j];
		}

	}

	return matrixHret;
}


#endif