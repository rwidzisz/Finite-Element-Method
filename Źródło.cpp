#include <iostream>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <vector>
#include <math.h>
#include <algorithm>

#include "Integral.h"
#include "Heat.h"
#include "Boundary.h"
#include "Pressure.h"

double H, W;
int nH, nW, conduct, ro, cp,alfa,aTemp,t0,deltaT,calT;
int global_integer = 1;

using namespace std;

void zczytaj_wymiary_siatki()
{
	int nr_linii = 1;
	string linia;

	fstream dane;
	dane.open("data.txt", ios::in);

	if (dane.good() == false)
	{
		cout << "Nie mozna otworzyc pliku";
		exit(0);
	}
	else
	{
		while (getline(dane, linia))
		{
			switch (nr_linii)
			{
			case 1:
				H = stof(linia);
				break;
			case 2:
				W = stof(linia);
				break;
			case 3:
				nH = stoi(linia);
				break;
			case 4:
				nW = stoi(linia);
				break;
			case 5:
				conduct = stoi(linia);
				break;
			case 6:
				ro = stoi(linia);
				break;
			case 7:
				cp = stoi(linia);
				break;
			case 8:
				alfa = stoi(linia);
				break;
			case 9:
				aTemp = stoi(linia);
				break;
			case 10:
				t0 = stoi(linia);
				break;
			case 11:
				deltaT = stoi(linia);
				break;
			case 12:
				calT = stoi(linia);
				break;
			}
			nr_linii++;
		}
	}
	dane.close();

}
 
std::vector<wezel_siatki> generuj_wezly()
{
	cout << endl << "Wezly" << endl;
	int nN = nH * nW;
	std::vector<wezel_siatki>zbior_wezlow;

	//Rozmiar wezlow
	double deltaX = (H / (nH - 1));
	double deltaY = (W / (nW - 1));

	//Inicjacja wezlow
	for (int i = 0; i < nW; i++) {
		for (int j = 0; j < nH; j++) {
			wezel_siatki *nowy = new wezel_siatki;
			nowy->x = (i*deltaX);
			nowy->y = (j*deltaY);
			if (nowy->x == 0 || i==nW-1 || nowy->y == 0 || j==nH-1 )
			{
				nowy->BC = true;
			}
			else {
				nowy->BC = false;
			}
			zbior_wezlow.push_back(*nowy);
		}
	}


	//Drukowanie
	for (int i = nH - 1; i >= 0; i--) {
		for (int j = 0; j < nW; j++) {
			cout << zbior_wezlow[i].x << "," << zbior_wezlow[j].y << "\t";
		}
		cout << endl << endl;
	}

	return zbior_wezlow;
}

std::vector<element_siatki> generuj_elementy(vector<wezel_siatki> x) {

	vector<wezel_siatki> zbior_wezlow = x;

	cout << endl << "Elementy siatki" << endl;
	int nE = (nH - 1) * (nW - 1);	//ilosc elementow
	std::vector<element_siatki>zbior_elementow;

	//Inicjacja elementow
	int element_integer = 0, check_column = 0, ID_integer = 1;
	for (int i = 0; i < nE; i++) {
		element_siatki *nowy = new element_siatki;

		if (check_column >= nH - 1)
		{
			element_integer++;
			ID_integer++;
			check_column = 0;
		}

		nowy->ID[0] = { ID_integer };
		nowy->ID[1] = { ID_integer + nH };
		nowy->ID[2] = { ID_integer + nH + 1 };
		nowy->ID[3] = { ID_integer + 1 };

		nowy->pkt1 = zbior_wezlow[element_integer];
		nowy->pkt2 = zbior_wezlow[element_integer + nH];
		nowy->pkt3 = zbior_wezlow[element_integer + nH + 1];
		nowy->pkt4 = zbior_wezlow[element_integer + 1];

		zbior_elementow.push_back(*nowy);

		element_integer++;
		ID_integer++;
		check_column++;
	}


	//Drukowanie elementow
	int nE_integer = 0;
	for (int i = nH - 1; i > 0; i--) {
		int ID_integer = i;
		for (int j = 0; j < nW - 1; j++) {
			cout << ID_integer << "," << ID_integer + nH << "," << ID_integer + nH + 1 << "," << ID_integer + 1;
			cout << " | ";
			ID_integer = ID_integer + nH;
		}
		cout << endl;
	}

	return zbior_elementow;
}


int main()
{
	integral integral;
	heat heat;
	boundary_condition bc;
	pressure pressure;

	//Siatka E 
	zczytaj_wymiary_siatki();
	
	int nt = calT / deltaT;
	int nE = (nH - 1) * (nW - 1);
	int nN = nH * nW;

	vector<wezel_siatki>zbior_wezlow;
	vector<element_siatki>zbior_elementow;

	zbior_wezlow = (generuj_wezly());
	zbior_elementow = (generuj_elementy(zbior_wezlow));

	//Punkty ca³kowania
	double p057 = sqrt(1.0 / 3.0);	//0.57
	double n057 = -1 * sqrt(1.0 / 3.0);	//-0.57
	pcalkowania dwu_pkt1(n057, n057), dwu_pkt2(p057, n057), dwu_pkt3(n057, p057), dwu_pkt4(p057, p057);

	double p077 = sqrt(3.0 / 5.0);	//0.77
	double n077 = -1 * p077;	//-0.77
	pcalkowania troj_pkt1(n077, n077), troj_pkt2(0, n077), troj_pkt3(p077, n077), troj_pkt4(n077, 0), troj_pkt5(0, 0), troj_pkt6(p077, 0), troj_pkt7(n077, p077), troj_pkt8(0, p077), troj_pkt9(p077, p077);

	double p086 = sqrt((3.0 / 7.0 + (2.0 / 7.0*sqrt(6.0 / 5.0))));
	double n086 = -1 * p086;
	double p033 = sqrt((3.0 / 7.0 - (2.0 / 7.0*sqrt(6.0 / 5.0))));
	double n033 = -1 * p033;
	pcalkowania czte_pkt1(n086, n086), czte_pkt2(n033, n086), czte_pkt3(p033, n086), czte_pkt4(p086, n086), czte_pkt5(n086, n033), czte_pkt6(n033, n033), czte_pkt7(p033, n033), czte_pkt8(p086, n033), czte_pkt9(n086, p033), czte_pkt10(n033, p033), czte_pkt11(p033, p033), czte_pkt12(p086, p033), czte_pkt13(n086, p086), czte_pkt14(n033, p086), czte_pkt15(p033, p086), czte_pkt16(p086, p086);

	//Tworzenie zbiorów punktów ca³kowania
	std::vector<pcalkowania>zbior_pcalkowania2;
	zbior_pcalkowania2.push_back(dwu_pkt1);
	zbior_pcalkowania2.push_back(dwu_pkt2);
	zbior_pcalkowania2.push_back(dwu_pkt3);
	zbior_pcalkowania2.push_back(dwu_pkt4);

	std::vector<pcalkowania>zbior_pcalkowania3;
	zbior_pcalkowania3.push_back(troj_pkt1);
	zbior_pcalkowania3.push_back(troj_pkt2);
	zbior_pcalkowania3.push_back(troj_pkt3);
	zbior_pcalkowania3.push_back(troj_pkt4);
	zbior_pcalkowania3.push_back(troj_pkt5);
	zbior_pcalkowania3.push_back(troj_pkt6);
	zbior_pcalkowania3.push_back(troj_pkt7);
	zbior_pcalkowania3.push_back(troj_pkt8);
	zbior_pcalkowania3.push_back(troj_pkt9);

	std::vector<pcalkowania>zbior_pcalkowania4;
	zbior_pcalkowania4.push_back(czte_pkt1);
	zbior_pcalkowania4.push_back(czte_pkt2);
	zbior_pcalkowania4.push_back(czte_pkt3);
	zbior_pcalkowania4.push_back(czte_pkt4);
	zbior_pcalkowania4.push_back(czte_pkt5);
	zbior_pcalkowania4.push_back(czte_pkt6);
	zbior_pcalkowania4.push_back(czte_pkt7);
	zbior_pcalkowania4.push_back(czte_pkt8);
	zbior_pcalkowania4.push_back(czte_pkt9);
	zbior_pcalkowania4.push_back(czte_pkt10);
	zbior_pcalkowania4.push_back(czte_pkt11);
	zbior_pcalkowania4.push_back(czte_pkt12);
	zbior_pcalkowania4.push_back(czte_pkt13);
	zbior_pcalkowania4.push_back(czte_pkt14);
	zbior_pcalkowania4.push_back(czte_pkt15);
	zbior_pcalkowania4.push_back(czte_pkt16);
	
	
	double**  matrixHG = new double*[nN];
	for (int i = 0; i < nN; ++i)
		matrixHG[i] = new double[nN] {0.0};

	double**  matrixCG = new double*[nN];
	for (int i = 0; i < nN; ++i)
		matrixCG[i] = new double[nN] {0.0};

	double**  matrixBC = new double*[nN];
	for (int i = 0; i < nN; ++i)
		matrixBC[i] = new double[nN] {0.0};

	double** vectorPC = new double*[nN];
	for (int i = 0; i < nN; ++i)
		vectorPC[i] = new double[nN] {0.0};
	
	double** matrixHprim = new double*[nN];
	for (int i = 0; i < nN; ++i)
		matrixHprim[i] = new double[nN] {0.0};

	double** vectorPprim = new double*[nN];
	for (int i = 0; i < nN; ++i)
		vectorPprim[i] = new double[nN] {0.0};

	static double vectorT[16][1];
	for(int i=0;i<16;i++)
	{
		vectorT[i][0] = t0;
	}

	int time = 50;
	
	//Zbiory przechowuj¹ce macierze i wektory lokalne
	std::vector<double**> zbior_mlokalnychH = std::vector<double**>();	
	std::vector<double**> zbior_mlokalnychC = std::vector<double**>();	
	std::vector<double**> zbior_mlokalnychBC = std::vector<double**>();	
	std::vector<double*> zbior_wlokalnychP = std::vector<double*>();
	

	cout << "\nRozwiazywanie ukladow rownan:" << endl;
	
	//Finalny wynik
	for (int i = 0; i < calT / deltaT; i++) {

		for (int number_element = 0; number_element < nE; number_element++) {	//Tworzenie macierzy
			zbior_mlokalnychH.push_back(integral.macierz_lokalna_dwupunkt(zbior_pcalkowania2, zbior_elementow[number_element], conduct));
			zbior_mlokalnychC.push_back(heat.macierz_lokalna_C_dwupunkt(zbior_pcalkowania2, zbior_elementow[number_element], cp, ro));
			zbior_mlokalnychBC.push_back(bc.macierz_lokalna_BC_dwupunkt(zbior_elementow[number_element], alfa));
			zbior_wlokalnychP.push_back(pressure.macierz_lokalna_P_dwupunkt(zbior_elementow[number_element], alfa, aTemp));

		}

		for (int number_element = 0; number_element < nE; number_element++) {	//Agregacja
			for (int i = 0; i < 4; i++) {
				for (int j = 0; j < 4; j++) {
					matrixHG[zbior_elementow[number_element].ID[i] - 1][zbior_elementow[number_element].ID[j] - 1] += zbior_mlokalnychH[number_element][i][j];
					matrixCG[zbior_elementow[number_element].ID[i] - 1][zbior_elementow[number_element].ID[j] - 1] += zbior_mlokalnychC[number_element][i][j];
					matrixBC[zbior_elementow[number_element].ID[i] - 1][zbior_elementow[number_element].ID[j] - 1] += zbior_mlokalnychBC[number_element][i][j];
				}
				vectorPC[0][zbior_elementow[number_element].ID[i] - 1] += zbior_wlokalnychP[number_element][i];
			}
		}

		//Uk³ady równañ

		//	HG+BC			
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				matrixHG[i][j] = matrixHG[i][j] + matrixBC[i][j];
			}
		}

		//	H+(C/dT)
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				matrixHprim[i][j] = matrixHG[i][j] + (matrixCG[i][j] / deltaT);
				//cout << matrixHpom[i][j] << " ";
			}
			//cout << endl;
		}

		cout << endl;
		//(CG/dT)
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				matrixCG[i][j] = (matrixCG[i][j] / deltaT);
			}
		}

		//CG*vectorT
		for (int i = 0; i < nN; ++i)		//rows1
			for (int j = 0; j < 1; ++j)		//colums2
				for (int k = 0; k < nN; ++k)//columns1
				{
					vectorPprim[i][j] += matrixCG[i][k] * vectorT[k][j];
				}

		//P-(C/dt)*vectorT
		for (int i = 0; i < nN; i++) {
			vectorPprim[i][0] = vectorPprim[i][0] - vectorPC[0][i];
			//cout << vectorPpom[i][0] << "	";
		}


		//Metoda Gauusa Jordana
		double matrixGJ[18][18];
		double x[18];
		double ratio;

		for (int i = 1; i <= 16; i++) {
			for (int j = 0; j <= 16 + 1; j++) {
				if (j == 17)
				{
					matrixGJ[i][j] = vectorPprim[i - 1][0];
				}
				else
				{
					matrixGJ[i][j] = matrixHprim[i - 1][j - 1];
				}
			}
		}

		for (int i = 1; i <= 16; i++)
		{
			if (matrixGJ[i][i] == 0.0)
			{
				cout << "Mathematical Error!";
				exit(0);
			}
			for (int j = 1; j <= 16; j++)
			{
				if (i != j)
				{
					ratio = matrixGJ[j][i] / matrixGJ[i][i];
					for (int k = 1; k <= 16 + 1; k++)
					{
						matrixGJ[j][k] = matrixGJ[j][k] - ratio * matrixGJ[i][k];
					}
				}
			}
		}

		/* Obtaining Solution */
		for (int i = 1; i <= 16; i++)
		{
			x[i] = matrixGJ[i][16 + 1] / matrixGJ[i][i];
		}
		
		//Przypisanie nowych wartoœæi do wektora T
		for (int i = 1; i <= 16; i++) {
			vectorT[i - 1][0] = x[i];
			//cout << vectorT[i-1][0] << "	";
		}

		double minT = vectorT[0][0];
		double maxT = vectorT[0][0];

		for (int i = 0; i < 16; i++) {
			if (minT > vectorT[i][0])
			{
				minT = vectorT[i][0];
			}
			if (maxT < vectorT[i][0])
			{
				maxT = vectorT[i][0];
			}
		}

		cout << endl << "Interacja:" << i << "	Time:" << time << "  MinT:" << minT << "  MaxT:" << maxT << endl;
		time += 50;

		//Zerowanie macierzy
		for (int i = 0; i < nN; i++) {
			for (int j = 0; j < nN; j++) {
				matrixHG[i][j] = NULL;
				matrixCG[i][j] = NULL;
				matrixBC[i][j] = NULL;
				vectorPC[i][j] = NULL;

				matrixHprim[i][j] = NULL;
				vectorPprim[j][0] = NULL;
			}
		}

		//Zerowanie zbiorów wektorów
		for (vector<double**>::iterator ne = zbior_mlokalnychH.begin(); ne != zbior_mlokalnychH.end(); ++ne) {
			delete **ne;
		}
		zbior_mlokalnychH.clear();

		for (vector<double**>::iterator ne = zbior_mlokalnychC.begin(); ne != zbior_mlokalnychC.end(); ++ne) {
			delete **ne;
		}
		zbior_mlokalnychC.clear();

		for (vector<double**>::iterator ne = zbior_mlokalnychBC.begin(); ne != zbior_mlokalnychBC.end(); ++ne) {
			delete **ne;
		}
		zbior_mlokalnychBC.clear();

		for (vector<double*>::iterator ne = zbior_wlokalnychP.begin(); ne != zbior_wlokalnychP.end(); ++ne) {
			delete *ne;
		}
		zbior_wlokalnychP.clear();

	}

	for (int i = 0; i < nN; ++i) {
		delete[] matrixHG[i];
	}
	delete[] matrixHG;

	for (int i = 0; i < nN; ++i) {
		delete[] matrixCG[i];
	}
	delete[] matrixCG;

	for (int i = 0; i < nN; ++i) {
		delete[] matrixBC[i];
	}
	delete[] matrixBC;

	for (int i = 0; i < 1; ++i) {
		delete[] vectorPC[i];

		delete[] vectorPC;
	}

	for (int i = 0; i < nN; ++i) {
		delete[] matrixHprim[i];
	}
	delete[] matrixHprim;


	for (int i = 0; i < 1; ++i) {
		delete[] vectorPprim[i];

		delete[] vectorPprim;
	}

	return 0;
}
