
#ifndef NACA4_H
#define NACA4_H

#include <vector>
#include <string>
#include "Airfoil.h"
#include "LinSpace.h"


/* Symbol definitions
 * a0l: zero-lift angle of attack
 * Cmac: prefer
 * c: chord
 * e: Oswald efficiency factor
 *
 * m: maximum thickness ordinate of the mean line in fraction of the chord
 *
 * p: x-location of maximum thickness, m
 * xu: abscissa points on the upper surface of the wing section
 * xl: abscissa points on the lower surface of the wing section
 * yu: ordinate points on the upper surface of the wing section
 * yl: ordinate points on the lower surface of the wing section
*/

class NACA4 : public Airfoil
{

public:
	NACA4(std::string NACA, int nseg = 128);

	double get_a0l();
	double getCmac();
	vector<double> getCamberLine();
	vector<double> getCamberLineSlope();

private:
	double p, m, e;
	vector<double> xc, yc, dycdx;
	vector<double> a0l, Cmac;

	int nseg;

	double calculate_a0l();
 	double calculateCmac();
	vector<double> calculateCamberLine();
	double calculateCamberLineSlope();
	double upperFunction(double abscissa);
	double lowerFunction(double ordinate);

	void createAbscissa();
	void createOrdinate();

};
#endif
