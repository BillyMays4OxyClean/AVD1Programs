#include <string>
using namespace std;

#ifndef AIRFOIL_H
#define AIRFOIL_H

/* Symbol definitions
 * a0l: zero-lift angle of attack
 * Cmac: prefer
 * c: chord
 * e: Oswald efficiency factor
 *
 * m: maximum thickness ordinate of the mean line in fraction of the chord
 *
 * p: x-location of maximum thickness chordwise position of m
 * xu: abscissa of point on the upper surface of a wing section
*/

class Airfoil
{
public:

	void Airfoil();

	double Calculate_a0l();

	double Calculate_Cmac();

	double CamberLine();

	double CamberLineSlope();

//	Airfoil ImportFromDatFile(const char * filename);

private:
	double m_p, m_m, m_t, m_c;
	double m_xu, m_xl, m_yu, m_yl, m_xc, m_yc, m_dycdx;
	double m_a0l, m_Cmac, m_e;

	string m_Name;

	int m_nseg;
};

#endif AIRFOIL_H
