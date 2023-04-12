#include <iostream>
#include <cmath>
#include <string>

#include "NACA4.h"

using namespace std;

NACA4::NACA4(){};

NACA4::NACA4(std::string NACA)
{
	/* Make a NACA string validator, based on a regex maybe. If it doesn't
	 pass the validator, an exception is made. */
	int p = NACA.at(0) - '0';
	int m = NACA.at(1) - '0';
	int t = stoi(NACA.substr(2,3));

      cout << "Generating NACA 4 digit airfoil: " << NACA << endl;
      cout << p << endl;
      cout << m << endl;
      cout << t << endl;
};

NACA4::NACA4(const char* NACA)
{
	/* Make a NACA string validator, based on a regex maybe. If it doesn't
	 pass the validator, an exception is made. */
	int p = NACA[0] - '0';
	int m = NACA[1] - '0';
	int t = stoi(string(NACA).substr(2,3));

      cout << "Generating NACA 4 digit airfoil: " << NACA << endl;
      cout << p << endl;
      cout << m << endl;
      cout << t << endl;
};