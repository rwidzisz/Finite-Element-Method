#pragma once
#include <vector>
#include "Mesh.h"
#include "Integral.h"

using std::vector;

class boundary_condition {
public:
	double **macierzHBC(double N1, double N2, double N3, double N4);
	double **macierz_lokalna_BC_dwupunkt(element_siatki element, int k);
	double **macierz_lokalna_BC_trzypunkt(element_siatki element, int k);
	double **macierz_lokalna_BC_czteropunkt(element_siatki element, int k);
};

