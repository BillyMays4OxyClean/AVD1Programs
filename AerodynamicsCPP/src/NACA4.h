#include "Airfoil.h"
#include <vector>
#include <string>

#ifndef NACA4_H
#define NACA4_H

class NACA4 : Airfoil
{

public:
	void NACA4();
	void NACA4(std::string NACA, double wingChord = 1, double oswaldFactor_e = 0.7);

	double get_a0l();
	double getCmac();
	vector<double> getCamberLine();
	vector<double> getCamberLineSlope();

private:
        double calculate_a0l();
        double calculateCmac();
        double calculateCamberLine();
        double calculateCamberLineSlope();

protected:

};
#endif
