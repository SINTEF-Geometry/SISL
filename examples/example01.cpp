/* SISL - SINTEF Spline Library version 4.4.                              */
/* Definition and interrogation of NURBS curves and surface.              */
/* Copyright (C) 1978-2005, SINTEF ICT, Applied Mathematics, Norway.      */

/* This program is free software; you can redistribute it and/or          */
/* modify it under the terms of the GNU General Public License            */
/* as published by the Free Software Foundation version 2 of the License. */

/* This program is distributed in the hope that it will be useful,        */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of         */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          */
/* GNU General Public License for more details.                           */

/* You should have received a copy of the GNU General Public License      */
/* along with this program; if not, write to the Free Software            */
/* Foundation, Inc.,                                                      */
/* 59 Temple Place - Suite 330,                                           */
/* Boston, MA  02111-1307, USA.                                           */

/* Contact information: e-mail: tor.dokken@sintef.no                      */
/* SINTEF ICT, Department of Applied Mathematics,                         */
/* P.O. Box 124 Blindern,                                                 */
/* 0314 Oslo, Norway.                                                     */

/* SISL commercial licenses are available for:                            */
/* - Building commercial software.                                        */
/* - Building software whose source code you wish to keep private.        */

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>

#include "sisl.h"
#include "GoReadWrite.h"

using namespace std;


namespace {
    string OUT_FILE_CURVE = "example1_curve.g2";
    string OUT_FILE_POINTS = "example1_points.g2";

    string DESCRIPTION = 
    "This program will generate a new SISL spline object from \n"
    "predefined control points and parametrization, and will \n"
    "write the result in Go-format to the file '" + OUT_FILE_CURVE + "'.\n"
    "The curve's control points will be written separately to \n"
    "the file '" + OUT_FILE_POINTS + "'.\n";

    const int number = 10;
    const int order = 4;
    
    double coef[] = {0,  0,   0,
		     1,  0, 0.5,
		     1,  1,   1,
		     0,  1, 1.5,
		     0,  0,   2,
		     1,  0, 2.5,
		     1,  1,   3,
		     0,  1, 3.5,
		     0,  0,   4,
		     1,  0, 4.5};

    double knots[] = {0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 7, 7, 7};

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {
	
	SISLCurve* curve = newCurve(number, // number of control points
				    order,  // order of spline curve (degree + 1)
				    knots,  // pointer to knot vector (parametrization)
				    coef,   // pointer to coefficient vector (control points)
				    1,      // kind = polynomial B-spline curve
				    3,      // dimension
				    0);     // no copying of information, 'borrow' arrays
	if (!curve) {
	    throw runtime_error("Error occured while generating curve.");
	}
	
	ofstream os_curve(OUT_FILE_CURVE.c_str());
	ofstream os_points(OUT_FILE_POINTS.c_str());
	if (!os_curve || !os_points) {
	    throw runtime_error("Unable to open output file.");
	}
	
	// write result to file
	writeGoCurve(curve, os_curve);
	writeGoPoints(number, coef, os_points);

	// cleaning up
	freeCurve(curve);
	os_curve.close();
	os_points.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
