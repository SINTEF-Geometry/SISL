/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

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
    // double epoint[12] = { 0, 0, 0,
    // 			  1, 0, 0,
    // 			  0, 0, 1,
    // 			  0.000000001, 1, 0 };
    // double eptyp[4] = { 1, 4, 1, 4 };
    s1606(cv[0], cv[1], epsge, point0, point1,
  	  blendtype, dim, order, &cv[2], &stat);
    // double astpar = 0.0;
    // double cendpar;
    // double aepsge = 1e-6;
    //    s1611(epoint, 4, 3, eptyp, 1, 4, astpar, aepsge, &cendpar, &cv[2], &stat);
    int n = cv[2]->in;
    int k = cv[2]->ik;
    // int t = cv[2]->ikind;
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
