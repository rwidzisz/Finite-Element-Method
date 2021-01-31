#pragma once
#include <vector>
#include "Mesh.h"
#include "Integral.h"

using std::vector;

class heat {
public:
	double *jakobian_C(double n, double E, element_siatki element);
	double **macierzC(double tab[5], int cp);
	double **macierz_lokalna_C_dwupunkt(vector<pcalkowania> ZPC, element_siatki element, int cp, int ro);
	double **macierz_lokalna_C_trzypunkt(vector<pcalkowania> ZPC, element_siatki element, int cp, int ro);
	double **macierz_lokalna_C_czteropunkt(vector<pcalkowania> ZPC, element_siatki element, int cp, int ro);
};