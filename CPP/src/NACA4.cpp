#include <iostream>
#include <cmath>

#include "NACA4.h"

using namespace std;

NACA4::NACA4(std::string NACA,int nseg=128)
:
c(1),
x(LinSpace(0,c,nseg))
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

vector<double>
NACA4::calculateCamberLine()
{
      /* Note to self: This c++ loop differs from its matlab counterpart in the sense that "i"
         i does not represent x in this case. It represents the current index of x.
         Therefore, another variable for x must be created. */
      
	for (size_t i = 0; i < (nseg-1); i++)
	{
            double xi = x.at(i);
            if (i <= p)
            {
                  yc.at(i) = m / pow(p,2) * (2*p*xi - pow(xi,2));
                  dycdx.at(i) = 2 * m / pow(p, 2.0) * (p - xi);
            }
            else if (xi > p)
            {
                  yc.at(i) = m / pow((1-p), 2.0) * ((1 - 2*p) + 2*p*xi - pow(xi,2));
                  dycdx.at(i) = 2 * m / pow((1-p), 2.0) * (p - xi);
            }
	}
}

double
NACA4::calculateCmac()
{
	
}
