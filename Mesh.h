#pragma once
#include <iostream>

class wezel_siatki {
public:
	double x;
	double y;
	bool BC;
};

class element_siatki {
public:
	wezel_siatki pkt1;
	wezel_siatki pkt2;
	wezel_siatki pkt3;
	wezel_siatki pkt4;

	int ID[4];
	element_siatki() {}

};
;