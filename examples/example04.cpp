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
    string IN_FILE_CURVE = "example1_curve.g2";
    string OUT_FILE_CURVE = "example4_curve.g2";

    string DESCRIPTION = 
    "This program will generate a SISL spline curve object that \n"
    "approximates an offset from a given input curve (note that \n"
    "the exact offset curve is not generally possible to \n"
    "represent exactly using splines).  The routine used is \n"
    "s1360.\n"
    "Input: " + IN_FILE_CURVE + " \n"
    "Output: " + OUT_FILE_CURVE + "\n\n";

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();


    try {
	ifstream is(IN_FILE_CURVE.c_str());
	if (!is) {
	    string error_message = 
		"Unable to open input file: " + IN_FILE_CURVE +
		".  Are you sure you have run the previous sample program?";
	    throw runtime_error(error_message.c_str());
	}

	SISLCurve* c1 = readGoCurve(is);

	double offset = double(1); // the offset between the 'old' curve and the result.
	double epsge = 1.0e-5;     // geometric tolerance 
	double direction[] = {0, 0, 1}; // the direction of the offset
	int dim = 3; // the dimension of the Euclidean space
	int jstat; // status variable
	SISLCurve* offset_curve = 0;

	s1360(c1,            // the 'old' curve
	      offset,        // the offset value
	      epsge,         // geometric tolerance
	      direction,     // offset direction
	      0,             // max step length.  0 indicate the longest box side of 's1'
	      dim,           // the dimension
	      &offset_curve, // the resulting offset curve
	      &jstat);       // status variable
	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1360.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1360. \n" << endl;
	}

	ofstream os(OUT_FILE_CURVE.c_str());
	if (!os) {
	    throw runtime_error("Unable to open output file.");
	}

	// write result to file
	writeGoCurve(offset_curve, os);
	
	// cleaning up
	freeCurve(c1);
	freeCurve(offset_curve);
	is.close();
	os.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
