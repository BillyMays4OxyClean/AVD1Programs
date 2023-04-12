#include "Airfoil.h"
#include <vector>
#include <string>

#ifndef NACA4_H
#define NACA4_H

class NACA4 : Airfoil
{

public:
	NACA4();
	NACA4(std::string NACA);
      NACA4(const char* NACA);

	double get_a0l();
	double getCmac();
	vector<double> getCamberLine();
	vector<double> getCamberLineSlope();

private:
      double calculate_a0l();
      double calculateCmac();
      double calculateCamberLine();
      double calculateCamberLineSlope();
      double upperFunction(double abscissa);
      double lowerFunction(double ordinate);


protected:

};
#endif
