/*
 * IntegrateSimpson.h
 *
 *  Created on: Nov 21, 2022
 *      Author: lcppa
 */

#ifndef INTEGRATESIMPSON_H_
#define INTEGRATESIMPSON_H_

/* Numerical Integration using Simpson's rule.*/

#include <cmath>
#include <iostream>
#include <vector>
using namespace std;

double IntegrateSimpson(vector<double> func, double a, double b)
{
   unsigned int nseg = func.size();

   if (pow(-1, nseg) < 0 || nseg < 2)
	   return printf("The input data array must be even to use Simpson's rule (nseg = %d).\n", nseg);

   double output{func.at(0)};

   for (unsigned int i{1}; i < nseg; i++)
   {
	   if (pow(-1, i) < 0 && i != nseg)
	   {
		   output += 4*func.at(i);
	   }
	   else if (pow(-1, i) > 0 && i != nseg && i != 0)
	   {
		   output += 2*func.at(i);
	   }
	   else if(i == nseg)
	   {
		   output += func.at(i);
	   };
   };

   double delta_x = (b - a)/nseg;
   output *= delta_x / 3;

   return output;
};

#endif /* INTEGRATESIMPSON_H_ */
