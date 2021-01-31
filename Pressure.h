#pragma once
#include <vector>
#include "Mesh.h"
#include "Integral.h"

using std::vector;

class pressure {
public:
	double *macierz_lokalna_P_dwupunkt(element_siatki element, int alfa, int talfa);
	double *macierz_lokalna_P_trójpunkt(element_siatki element, int alfa, int talfa);
	double *macierz_lokalna_P_czteropunkt(element_siatki element, int alfa, int talfa);
};
