//===========================================================================
//                                                                           
// File: s1606prob.C                                                         
//                                                                           
// Created: Wed Aug 30 13:46:14 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: s1606prob.C,v 1.1 2001-09-05 14:30:12 afr Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#include "sisl.h"
#include <iostream>


using namespace std;


int main()
{
    double knots[4]
	= {0.0, 0.0, 1.0, 1.0};
    double c0[6]
	= { 0.000002, 0.000150, 0.000000, 0.000000, 0.000000, 0.000000 };
    double c1[6]
	= { 0.000000, 0.001768, -0.001768, 0.000000, 0.001914, -0.001914 };
//      double c0[6]
//  	= { 0, 0, 0, 1, 0, 0 };
//      double c1[6]
//  	= { 2, 1, 1, 2, 2, 1 };

    SISLCurve* cv[3];
    cv[0] = newCurve(2, 2, knots, c0, 1, 3, 0);
    cv[1] = newCurve(2, 2, knots, c1, 1, 3, 0);
    cv[2] = 0;
    double epsge = 1e-6;
    double point0[3] = {0.000000, 0.001768, -0.001768};
    double point1[3] = {0.000002, 0.000150,  0.000000};
//      double point0[3] = { 1, 0, 0 };
//      double point1[3] = { 2, 1, 1 };
    int blendtype = 2;
    int dim = 3;
    int order = 4;
    int stat;
//      double epoint[12] = { 0.000002, 0.000150, 0.000000,
//  			  0.013332148306149432, 0.99991112296120743, 0,
//  			  0, 0.001768, -0.001768,
//  			  0, 0.70710678118654746, -0.70710678118654746 };
    double epoint[12] = { 0, 0, 0,
			  1, 0, 0,
			  0, 0, 1,
			  0.000000001, 1, 0 };
    double eptyp[4] = { 1, 4, 1, 4 };
    s1606(cv[0], cv[1], epsge, point0, point1,
  	  blendtype, dim, order, &cv[2], &stat);
    double astpar = 0.0;
    double cendpar;
    double aepsge = 1e-6;
    //    s1611(epoint, 4, 3, eptyp, 1, 4, astpar, aepsge, &cendpar, &cv[2], &stat);
    int n = cv[2]->in;
    int k = cv[2]->ik;
    int t = cv[2]->ikind;
    int d = cv[2]->idim;
    cout.precision(15);
    cout << "GoNurbsCurve3D\ndimension 3\nbounding_box 0\n"
	 << "attribute_list_length 0\nparameters "
	 << cv[2]->et[k-1] << ' ' << cv[2]->et[n] << '\n'
	 << "periodic 0\nbspline_basis\nnum_coeffs "
	 << n << "\norder " << k << "\nknot_vector\n";
    for (int i = 0; i < n + k; ++i) {
  	cout << cv[2]->et[i] << ' ';
    }
    cout << "\nkind 1\ndim 3\ncoefficients\n";
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < d; ++j) {
	    cout << cv[2]->ecoef[i*d + j] << ' ';
	}
	cout << '\n';
    }
    cout << "closedness -1" << endl;
}
