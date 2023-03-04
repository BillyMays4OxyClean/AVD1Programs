#include "NACA4.h"
#include <cmath>
#include <string>
void
NACA4::NACA4
(std::string NACA, double wingChord = 1.0, double oswaldFactor_e = 0.7)
{
	/* Make a NACA string validator, based on a regex maybe. If it doesn't
	 pass the validator, an exception is made. */
	int p = std::stoi (static_cast<std::string>(NACA.at(0)));
	int m = std::stoi (static_cast<std::string>(NACA.at(0)));
	int t = std::stoi (static_cast<std::string>(NACA.at(0)));
};
