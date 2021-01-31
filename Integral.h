#pragma once
#include <vector>
#include "Mesh.h"
#include "Integral.h"

using std::vector;

class pcalkowania
{
public:
	double x;
	double y;

	pcalkowania(double xx, double yy)
	{
		this->x = xx;
		this->y = yy;
	}
};

class integral {
public:
	double *jakobian(double n, double E, element_siatki element);
	double **macierzH(double tab[9],int k);
	double **macierz_lokalna_dwupunkt(vector<pcalkowania> ZPC, element_siatki element, int k);
	double **macierz_lokalna_trzypunkt(vector<pcalkowania> ZPC, element_siatki element, int k);
	double **macierz_lokalna_czteropunkt(vector<pcalkowania> ZPC, element_siatki element, int k);

};