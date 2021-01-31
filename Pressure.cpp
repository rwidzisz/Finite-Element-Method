#ifndef PRESSURE_H
#define PRESSURE_H

#include "Pressure.h"

using namespace std;

double *pressure::macierz_lokalna_P_dwupunkt(element_siatki element, int alfa, int talfa) {

	double p057 = sqrt(1.0 / 3.0);	//0.57
	double n057 = -1 * p057;
	pcalkowania p1(n057, -1), p2(p057, -1), p3(1, n057), p4(1, p057), p5(p057, 1), p6(n057, 1), p7(-1, p057), p8(-1, n057);


	double matrixP[1][4];
	for (int i = 0; i < 1; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixP[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N1p1 = 0.5*(1 - p1.x);
		double N2p1 = 0.5*(1 + p1.x);
	
		double N1p2 = 0.5*(1 - p2.x);
		double N2p2 = 0.5*(1 + p2.x);

		double vector1[1][4] = { N1p1, N2p1,0,0 };
		double vector2[1][4] = { N1p2, N2p2,0,0 };

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] + vector2[0][j])*talfa*det);
		}


	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N2p3 = 0.5*(1 - p3.y);
		double N3p3 = 0.5*(1 + p3.y);

		double N2p4 = 0.5*(1 - p4.y);
		double N3p4 = 0.5*(1 + p4.y);

		double vector1[1][4] = { 0,N2p3,N3p3,0 };
		double vector2[1][4] = { 0,N2p4,N3p4,0 };

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] + vector2[0][j])*talfa*det);
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{
		double N3p5 = 0.5*(1 - p5.x);
		double N4p5 = 0.5*(1 + p5.x);
		 
		double N3p6 = 0.5*(1 - p6.x);
		double N4p6 = 0.5*(1 + p6.x);

		double vector1[1][4] = { 0,0,N3p5,N4p5 };
		double vector2[1][4] = { 0,0,N3p6,N4p6 };

		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] + vector2[0][j])*talfa*det);
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true)
	{
		double N1p7 = 0.5*(1 + p7.y);
		double N4p7 = 0.5*(1 - p7.y);

		double N1p8 = 0.5*(1 + p8.y);
		double N4p8 = 0.5*(1 - p8.y);

		double vector1[1][4] = { N1p7,0,0,N4p7 };
		double vector2[1][4] = { N1p8,0,0,N4p8 };

		double det = (element.pkt4.y - element.pkt1.y) / 2.0;
		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1*alfa*(vector1[0][j] + vector2[0][j])*talfa*det);		
		}
	}

	/*cout << "\n\nLokalna Macierz P" << endl;
	for (int i = 0; i < 1; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixP[i][j];

			if (j == 4 - 1)
				cout << endl;
		}*/


	double *matrixPret = new double[4];
	matrixPret[0] = matrixP[0][0];
	matrixPret[1] = matrixP[0][1];
	matrixPret[2] = matrixP[0][2];
	matrixPret[3] = matrixP[0][3];

	return matrixPret;
}

double *pressure::macierz_lokalna_P_trójpunkt(element_siatki element, int alfa, int talfa) {

	double p077 = sqrt(3.0 / 5.0);	//0.77
	double n077 = -1 * p077;	//-0.77
	pcalkowania p1(n077, -1), p2(0, -1), p3(p077, -1),
				p4(-1, n077), p5(-1, 0), p6(-1, p077),
				p7(p077, 1), p8(0, 1), p9(n077, 1),
				p10(-1, p077), p11(-1, 0), p12(-1, n077);

	double wg[3] = { (5.0 / 9.0),	(8.0 / 9.0), (5.0 / 9.0) };

	double matrixP[1][4];
	for (int i = 0; i < 1; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixP[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N1p1 = 0.5*(1 - p1.x);
		double N2p1 = 0.5*(1 + p1.x);

		double N1p2 = 0.5*(1 - p2.x);
		double N2p2 = 0.5*(1 + p2.x);

		double N1p3 = 0.5*(1 - p3.x);
		double N2p3 = 0.5*(1 + p3.x);

		double vector1[1][4] = { N1p1, N2p1,0,0 };
		double vector2[1][4] = { N1p2, N2p2,0,0 };
		double vector3[1][4] = { N1p3, N2p3,0,0 };

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j]*wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2])*talfa*det);
		}


	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N2p4 = 0.5*(1 - p4.y);
		double N3p4 = 0.5*(1 + p4.y);

		double N2p5 = 0.5*(1 - p5.y);
		double N3p5 = 0.5*(1 + p5.y);

		double N2p6 = 0.5*(1 - p6.y);
		double N3p6 = 0.5*(1 + p6.y);

		double vector1[1][4] = { 0,N2p4,N3p4,0 };
		double vector2[1][4] = { 0,N2p5,N3p5,0 };
		double vector3[1][4] = { 0,N2p6,N3p6,0 };

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2])*talfa*det);
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{
		double N3p7 = 0.5*(1 - p7.x);
		double N4p7 = 0.5*(1 + p7.x);

		double N3p8 = 0.5*(1 - p8.x);
		double N4p8 = 0.5*(1 + p8.x);

		double N3p9 = 0.5*(1 - p9.x);
		double N4p9 = 0.5*(1 + p9.x);

		double vector1[1][4] = { 0,0,N3p7,N4p7 };
		double vector2[1][4] = { 0,0,N3p8,N4p8 };
		double vector3[1][4] = { 0,0,N3p9,N4p9 };
	

		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{

			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2])*talfa*det);
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true)
	{
		double N1p10 = 0.5*(1 + p10.y);
		double N4p10 = 0.5*(1 - p10.y);

		double N1p11 = 0.5*(1 + p11.y);
		double N4p11 = 0.5*(1 - p11.y);

		double N1p12 = 0.5*(1 + p12.y);
		double N4p12 = 0.5*(1 - p12.y);

		double vector1[1][4] = { N1p10,0,0,N4p10 };
		double vector2[1][4] = { N1p11,0,0,N4p11 };
		double vector3[1][4] = { N1p12,0,0,N4p12 };

		double det = (element.pkt4.y - element.pkt1.y) / 2.0;
		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2])*talfa*det);
		}
	}

	/*cout << "\n\nLokalna Macierz P" << endl;
	for (int i = 0; i < 1; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixP[i][j];

			if (j == 4 - 1)
				cout << endl;
		}*/


	double *matrixPret = new double[4];
	matrixPret[0] = matrixP[0][0];
	matrixPret[1] = matrixP[0][1];
	matrixPret[2] = matrixP[0][2];
	matrixPret[3] = matrixP[0][3];

	return matrixPret;
}

double *pressure::macierz_lokalna_P_czteropunkt(element_siatki element, int alfa, int talfa) {

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

	double matrixP[1][4];
	for (int i = 0; i < 1; ++i) {
		for (int j = 0; j < 4; ++j)
		{
			matrixP[i][j] = 0;
		}
	}

	if (element.pkt1.BC == true && element.pkt2.BC == true)
	{
		double N1p1 = 0.5*(1 - p1.x);
		double N2p1 = 0.5*(1 + p1.x);

		double N1p2 = 0.5*(1 - p2.x);
		double N2p2 = 0.5*(1 + p2.x);

		double N1p3 = 0.5*(1 - p3.x);
		double N2p3 = 0.5*(1 + p3.x);

		double N1p4 = 0.5*(1 - p4.x);
		double N2p4 = 0.5*(1 + p4.x);

		double vector1[1][4] = { N1p1, N2p1,0,0 };
		double vector2[1][4] = { N1p2, N2p2,0,0 };
		double vector3[1][4] = { N1p3, N2p3,0,0 };
		double vector4[1][4] = { N1p4, N2p4,0,0 };

		double det = (element.pkt2.x - element.pkt1.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2]+vector4[0][j] * wg[3])*talfa*det);
		}


	}
	if (element.pkt2.BC == true && element.pkt3.BC == true)
	{
		double N2p5 = 0.5*(1 - p5.y);
		double N3p5 = 0.5*(1 + p5.y);

		double N2p6 = 0.5*(1 - p6.y);
		double N3p6 = 0.5*(1 + p6.y);

		double N2p7 = 0.5*(1 - p7.y);
		double N3p7 = 0.5*(1 + p7.y);

		double N2p8 = 0.5*(1 - p8.y);
		double N3p8 = 0.5*(1 + p8.y);

		double vector1[1][4] = { 0,N2p5,N3p5,0 };
		double vector2[1][4] = { 0,N2p6,N3p6,0 };
		double vector3[1][4] = { 0,N2p7,N3p7,0 };
		double vector4[1][4] = { 0,N2p8,N3p8,0 };

		double det = (element.pkt3.y - element.pkt2.y) / 2.0;

		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2] + vector4[0][j] * wg[3])*talfa*det);
		}
	}
	if (element.pkt3.BC == true && element.pkt4.BC == true)
	{
		double N3p9 = 0.5*(1 - p9.x);
		double N4p9 = 0.5*(1 + p9.x);
		
		double N3p10 = 0.5*(1 - p10.x);
		double N4p10 = 0.5*(1 + p10.x);
		
		double N3p11 = 0.5*(1 - p11.x);
		double N4p11 = 0.5*(1 + p11.x);
		
		double N3p12 = 0.5*(1 - p12.x);
		double N4p12 = 0.5*(1 + p12.x);


		double vector1[1][4] = { 0,0,N3p9,N4p9 };
		double vector2[1][4] = { 0,0,N3p10,N4p10 };
		double vector3[1][4] = { 0,0,N3p11,N4p11 };
		double vector4[1][4] = { 0,0,N3p12,N4p12 };
	


		double det = (element.pkt3.x - element.pkt4.x) / 2.0;

		for (int j = 0; j < 4; ++j)
		{

			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2] + vector4[0][j] * wg[3])*talfa*det);
		}
	}
	if (element.pkt4.BC == true && element.pkt1.BC == true)
	{
		double N1p13 = 0.5*(1 + p13.y);
		double N4p13 = 0.5*(1 - p13.y);

		double N1p14 = 0.5*(1 + p14.y);
		double N4p14 = 0.5*(1 - p14.y);

		double N1p15 = 0.5*(1 + p15.y);
		double N4p15 = 0.5*(1 - p15.y);

		double N1p16 = 0.5*(1 + p16.y);
		double N4p16 = 0.5*(1 - p16.y);

		double vector1[1][4] = { N1p13,0,0,N4p13 };
		double vector2[1][4] = { N1p14,0,0,N4p14 };
		double vector3[1][4] = { N1p16,0,0,N4p16 };
		double vector4[1][4] = { N1p16,0,0,N4p16 };
		
		double det = (element.pkt4.y - element.pkt1.y) / 2.0;
		for (int j = 0; j < 4; ++j)
		{
			matrixP[0][j] = matrixP[0][j] + (-1 * alfa*(vector1[0][j] * wg[0] + vector2[0][j] * wg[1] + vector3[0][j] * wg[2] + vector4[0][j] * wg[3])*talfa*det);
		}
	}

	/*cout << "\n\nLokalna Macierz P" << endl;
	for (int i = 0; i < 1; ++i)
		for (int j = 0; j < 4; ++j)
		{
			cout << " " << matrixP[i][j];

			if (j == 4 - 1)
				cout << endl;
		}*/


	double *matrixPret = new double[4];
	matrixPret[0] = matrixP[0][0];
	matrixPret[1] = matrixP[0][1];
	matrixPret[2] = matrixP[0][2];
	matrixPret[3] = matrixP[0][3];

	return matrixPret;
}



#endif