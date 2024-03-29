#include <iostream>
#include <cmath>

#include "NACA4.h"

using namespace std;

NACA4::NACA4(std::string NACA,int nseg):
xc(LinSpace(0.0,c,nseg)),
dycdx(LinSpace(0.0,0.0,nseg))
{
	/* Make a NACA string validator, based on a regex maybe. If it doesn't
	 pass the validator, an exception is made. */
	
	Name = NACA;
	p = static_cast<double>(NACA.at(0) - '0') / 10;
	m = static_cast<double>(NACA.at(1) - '0') / 100;
	t = stod(NACA.substr(2,3)) / 100;

	cout << "Generating NACA4 digit airfoil: " << NACA << endl;

      x = LinSpace(0.0,c,nseg);
      yu = LinSpace(0.0,c,nseg);
      yl = LinSpace(0.0,c,nseg);
      yc = LinSpace(0.0,c,nseg);
	vector<double> yt{LinSpace(0.0,0.0,nseg)};
	vector<double> theta{LinSpace(0.0,0.0,nseg)};

	for (size_t i = 0; i < nseg ; i++)
	{
		double xi = x.at(i);
		if (xi <= p)
		{
			yc.at(i) = m / pow(p,2) * (2*p*xi - pow(xi,2));
			dycdx.at(i) = 2 * m / pow(p, 2.0) * (p - xi);
		}
		else if (xi >= p)
		{
			yc.at(i) = m / pow((1-p), 2.0) * ((1 - 2*p) + 2*p*xi - pow(xi,2.0));
			dycdx.at(i) = 2 * m / pow((1-p), 2.0) * (p - xi);
		}

		theta.at(i) = atan(dycdx.at(i));

		yt.at(i) = t/0.2 * ( 0.29690*sqrt(xi) - 0.12600*xi - 0.35160 * pow(xi,2) + 0.28430 * pow(xi,3) - 0.10150 * pow(xi,4));

		if ((p != 0) || (m != 0))
		{
			yu.at(i) = yc.at(i) + yt.at(i)*cos(theta.at(i));
			yl.at(i) = yc.at(i) - yt.at(i)*cos(theta.at(i));
		}
		else
		{
			yu.at(i) = yt.at(i);
			yl.at(i) = -yt.at(i);
		}
		cout << yl.at(i) << ", " << yu.at(i) << endl;
	}

	ExportFile("NACA" + Name + ".dat");
};
