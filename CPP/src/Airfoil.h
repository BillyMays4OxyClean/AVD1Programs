#include <string>
#include <vector>
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
 * p: x-location of maximum thickness, m
 * xu: abscissa points on the upper surface of the wing section
 * xl: abscissa points on the lower surface of the wing section
 * yu: ordinate points on the upper surface of the wing section
 * yl: ordinate points on the lower surface of the wing section
*/

class Airfoil
{
public:

	void Airfoil();

	int ExportDatFile(std::string filename);
	int ImportFromDatFile(std::string filename);

	double getMaxThicknessLocation();
	double getMaxThickness();
	double get

	void PlotAirfoil();

private:
	double m_p, m_m, m_t, m_c;
	vector<double> m_xu, m_xl, m_yu, m_yl, m_xc, m_yc, m_dycdx;
	vector<double> m_a0l, m_Cmac, m_e;

	string m_Name;

	int m_nseg;
};

#endif AIRFOIL_H
