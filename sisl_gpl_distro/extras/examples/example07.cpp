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
    string IN_FILE_CURVE = "example6_curve_1.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program computes the length of a curve read from file.\n"
    "The name of this curve is: '" + IN_FILE_CURVE + "'.\n\n";

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
	    throw runtime_error("Unable to open input file.");
	}

	// read curve from file
	SISLCurve* cv = readGoCurve(is);

	// computing length of curve
	double epsge = 1.0e-5; // geometric tolerance
	double length = 0;
	int jstat = 0;

	s1240(cv,      // the curve we want to know the length of
	      epsge,   // geometric tolerance
	      &length, // the calculated length
	      &jstat); // status report variable

	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1240.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1240. \n" 
		 << endl;
	}

	cout << "Computed length of curve: " << length << "\n";

	// cleaning up
	freeCurve(cv);
	is.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
