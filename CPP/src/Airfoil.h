#include <string>
#include <vector>
using namespace std;

#include "LinSpace.h"

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
	Airfoil();
	
	int ExportFile(std::string filename);
      int ExportFile(std::string filename, std::vector<double> abscissa, std::vector<double> ordinate);
	int ImportFromDatFile(std::string filename);

	double getMaxThicknessLocation();
	double getMaxThickness();

	void PlotAirfoil();

	std::string Name;
	
	vector<double> x, yu, yl, yc;

	int nseg{128};
	double t, c{1.0};

};

#endif
