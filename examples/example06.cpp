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
    string OUT_FILE_CURVE_1 = "example6_curve_1.g2";
    string OUT_FILE_CURVE_2 = "example6_curve_2.g2";
    string OUT_FILE_POINT   = "example6_isectpoints.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program generates two curves (from internal data), and \n"
    "computes their intersections, using SISL routine s1857.  \n"
    "The curves are planar, but lying in 3D space.  The curves \n"
    "and the intersections will then be saved to file.  The \n"
    "curves will be saved to the two files '" + OUT_FILE_CURVE_1 
    + "' \nand '" + OUT_FILE_CURVE_2 + "'. The intersection points \n"
    " will be saved to the file '" + OUT_FILE_POINT + "'\n\n";

    const int dim = 3;

    const int c1_number = 4;
    const int c1_order = 4;
    
    double c1_coef[] = {0, 0, 0,
			0, 2, 0,
			2, 2, 0,
			2, 0, 0};
    double c1_knots[] = {0, 0, 0, 0, 1, 1, 1, 1};

    const int c2_number = 5;
    const int c2_order = 3;
    
    double c2_coef[] = {  0,   1,  0,
			0.5, 0.5,  0,
			  1,   2,  0,
			1.5, 0.5,  0,
			  2,   1,  0};
			
    double c2_knots[] = {0, 0, 0, 1, 2, 3, 3, 3};

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {
	ofstream os_cv1(OUT_FILE_CURVE_1.c_str());
	ofstream os_cv2(OUT_FILE_CURVE_2.c_str());
	ofstream os_pts(OUT_FILE_POINT.c_str());
	if (!os_cv1 || !os_cv2 || !os_pts) {
	    throw runtime_error("Unable to open output file.");
	}

	// generating curves from internal data
	SISLCurve* c1 = newCurve(c1_number, // num. of control points
				 c1_order,  // spline order
				 c1_knots,  // knotvector
				 c1_coef,   // control points
				 1,         // kind = polynomial B-spline curve
				 dim,       // dimension of space (3D)
				 1);        // copy input arrays

	SISLCurve* c2 = newCurve(c2_number, // num. of control points
				 c2_order,  // spline order
				 c2_knots,  // knotvector
				 c2_coef,   // control points
				 1,         // kind = polynomial B-spline curve
				 dim,       // dimension of space (3D)
				 1);        // copy input arrays
		 
	if(!c1 || !c2) {
	    throw runtime_error("Error occured while generating curves.");
	}

	// calculating intersection points
	double epsco = 1.0e-15; // computational epsilon
	double epsge = 1.0e-5; // geometric tolerance
	int num_int_points = 0; // number of found intersection points
	double* intpar1 = 0; // parameter values for the first curve in the intersections
	double* intpar2 = 0; // parameter values for the second curve in the intersections
	int num_int_curves = 0;   // number of intersection curves
	SISLIntcurve** intcurve = 0; // pointer to array of detected intersection curves
	int jstat; // status variable

	s1857(c1,              // first curve 
	      c2,              // second curve
	      epsco,           // computational resolution
	      epsge,           // geometry resolution
	      &num_int_points, // number of single intersection points
	      &intpar1,        // pointer to array of parameter values
	      &intpar2,        //               "
	      &num_int_curves, // number of detected intersection curves
	      &intcurve,       // pointer to array of detected intersection curves.
	      &jstat);
	
	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1857.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1857. \n" 
		     << endl;
	}

	// In this example, we do not expect to find intersection curves, since the
	// input curves do not overlap for more than one point at a time.

	cout << "Number of intersection points detected: " << num_int_points << endl;
	cout << "Number of intersection curves detected: " << num_int_curves << endl;
	
	// evaluating intersection points and writing them to file
	vector<double> point_coords_3D(3 * num_int_points);
	int i;
	for (i = 0; i < num_int_points; ++i) {
	    // calculating position, using curve 1
	    // (we could also have used curve 2, which would give approximately
	    // the same points).
	    int temp;
	    s1227(c1,         // we evaluate on the first curve
		  0,          // calculate no derivatives
		  intpar1[i], // parameter value on which to evaluate
		  &temp,      // not used for our purposes (gives parameter interval)
		  &point_coords_3D[3 * i], // result written here
		  &jstat);

	    if (jstat < 0) {
		throw runtime_error("Error occured inside call to SISL routine s1227.");
	    } else if (jstat > 0) {
		cerr << "WARNING: warning occured inside call to SISL routine s1227. \n" 
		     << endl;
	    }
	}

	// writing curves and intersections to file
	writeGoCurve(c1, os_cv1);
	writeGoCurve(c2, os_cv2);
	writeGoPoints(num_int_points, &point_coords_3D[0], os_pts);

	// cleaning up
	freeCurve(c1);
	freeCurve(c2);
	os_cv1.close();
	os_cv2.close();
	os_pts.close();
	free(intpar1);
	free(intpar2);
	for (i = 0; i < num_int_curves; ++i) {
	    freeIntcurve(intcurve[i]);
	}
	free(intcurve);

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
