
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

class NACA4 : Airfoil
{

public:
	NACA4(std::string NACA,int nseg=128);
	NACA4(const char* NACA,int nseg=128);

	double get_a0l();
	double getCmac();
	vector<double> getCamberLine();
	vector<double> getCamberLineSlope();

private:
	double m_p, m_m, m_t, m_c, m_e;
	vector<double> m_x, m_yu, m_yl, m_xc, m_yc, m_dycdx;
	vector<double> m_a0l, m_Cmac, m_e;

     	double calculate_a0l();
 	double calculateCmac();
        vector<double> calculateCamberLine();
        double calculateCamberLineSlope();
        double upperFunction(double abscissa);
        double lowerFunction(double ordinate);

	void createAbscissa();
	void createOrdinate();


protected:

};
#endif
