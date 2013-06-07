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

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstdlib>

#include "sisl.h"
#include "GoReadWrite.h"

using namespace std;

namespace {
    string OUT_FILE_CURVES = "example10_curves.g2";
    string OUT_FILE_SURFACE = "example10_surf.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program generates a series of curves, and then  \n"
    "generates a lofted surface through these curves.  The \n"
    "routine used is s1538.  The resulting curves and surfaces \n"
    " are written respectively to '" + OUT_FILE_CURVES + "' and\n"
    "'" + OUT_FILE_SURFACE + "'.\n\n";
    //==========================================================

    const int num_curves = 16;
    const int num_control_points = 4;
    const int curve_order = 4;
    const int loft_order = 4;
    
    // defining coefficients for the first of the curves
    double curve_coef[] = {0, 0, 0,
			   1, 2, 0,
			   2, 3, 0,
			   3, 0, 0};

    // defining the knotvector of the curve
    double curve_kvec[] = {0, 0, 0, 0, 1, 1, 1, 1};

    // values used to successfully change the curve shape
    const double x_dilate = 1.15;
    const double y_increment = 2;
    const double z_increment = 2;
    
    void incrementally_change_control_points();

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {
	ofstream os_curves(OUT_FILE_CURVES.c_str());
	ofstream os_surf(OUT_FILE_SURFACE.c_str());
	if (!os_surf || !os_curves) {
	    throw runtime_error("Unable to open output file.");
	}

	// generating series of curves
	vector<SISLCurve*> curves(num_curves);
	for (int c = 0; c < num_curves; ++c) {
	    curves[c] = newCurve(num_control_points, // number of curve control points
				 curve_order,        // the order of the curve 
				 curve_kvec,         // pointer to the knotvector
				 curve_coef,         // pointer to the control points,
				 1,                  // kind = 1: polynomial spline curve
				 3,                  // dimension of Euclidean space
				 1);                 // copy input arrays
	    if (!curves[c]) {
		throw runtime_error("Unable to generate curve.");
	    }

	    incrementally_change_control_points();

	    // writing curve to file
	    writeGoCurve(curves[c], os_curves);
	}

	// generating lofted surface
	SISLSurf* result_surf = 0;
	double* gpar = 0;
	vector<int> cv_type(num_curves, 1); // to indicate that all curves are 'ordinary'
	int jstat;
	
	s1538(num_curves,   // the number of curves that are used to generate the lofted surf.
	      &curves[0],   // pointer to the array of input curves
	      &cv_type[0],  // indicate the type of the input curves
	      0,            // start parameter for lofting direction
	      1,            // flag telling that the surface should be open
	      loft_order,   // order of the surface in the lofting direction
	      0,            // do not adjust tangents in the derivative curve
	      &result_surf, // the resulting surface
	      &gpar,        // get the parametrization along the lofted dir 
	      &jstat);      // status variable
	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1538.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1538. \n" 
		 << endl;
	}

	// writing surface to file
	writeGoSurface(result_surf, os_surf);

	// cleaning up
	free(gpar);
	for (int i = 0; i < num_curves; ++i) {
	    freeCurve(curves[i]);
	}
	freeSurf(result_surf);
	os_surf.close();
	os_curves.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};

namespace {

void incrementally_change_control_points()
{
    // incremental change in control points, so that the next curve will differ
    // (Feel free to experiment with your own formulas here!)

    for (int p = 0; p < num_control_points; ++p) {
	// incrementing z value
	curve_coef[3 * p + 0] *= x_dilate;
	curve_coef[3 * p + 2] += z_increment;
	    }
    curve_coef[3 * 1 + 1] -= 0.1 * y_increment;  // subtract increment value here
    curve_coef[3 * 2 + 1] += y_increment;  // add the increment value here
}

}; // end anonymous namespace
