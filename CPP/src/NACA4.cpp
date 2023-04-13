#include <iostream>
#include <cmath>
#include <string>

#include "NACA4.h"

using namespace std;

NACA4::NACA4(std::string NACA,int nseg=128)
:
c(1),
m_x(LinSpace(0,c,nseg))
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

	createAbscissa();	

};

NACA4::NACA4(const char* NACA,int nseg=128)
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

void
NACA4::createAbscissa()
{
	m_x = LinSpace(0,c,nseg);
}

vector<double>
NACA4::calculateCamberLine()
{
	for (i = 0; i < (nseg-1); i++)
	{

	}
}

double
NACA4::calculateCmac()
{
	
}
