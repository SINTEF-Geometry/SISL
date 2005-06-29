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
    string OUT_FILE_CURVE = "example5_curve.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program generates some conic sections represented as \n"
    "SISLCurves.  Since the NURBS representation is used, the \n"
    "generated curves would be exact. The SISL routine used is\n"
    "s1011.  The resulting curves will be saved in Go-format to \n"
    "the file " + OUT_FILE_CURVE + ".\n\n";

    const int num_shapes = 7; // number of generated shapes 
    const int dim = 3;

    // logically, the following arrays are 'const'.  But this generates calling error
    // with SISL,so we leave them un-consted for now.
    double shape[] = {0.01,  // an ellipse
		      0.25,  // an ellipse
		      0.4,   // an ellipse
		      0.5,   // a parabola
		      0.6,   // a hyperbola
		      0.75,  // a hyperbola
		      0.99}; // a hyperbola
    double start_pos[] = {0, 0, 0}; // start position (3D)
    double top_pos[] = {1, 2, 0};   // shoulder point (3D)
    double end_pos[] = {2, 0, 0};   // end position

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {

	ofstream os(OUT_FILE_CURVE.c_str());
	if (!os) {
	    throw runtime_error("Unable to open output file.");
	}

	for (int i = 0; i < num_shapes; ++i) {

	    SISLCurve* result_curve = 0;
	    int jstat;
	    
	    // generate the conic section
	    s1011(start_pos,
		  top_pos,
		  end_pos,
		  shape[i],
		  dim,
		  &result_curve,
		  &jstat);
	    if (jstat < 0) {
		throw runtime_error("Error occured inside call to SISL routine s1011.");
	    } else if (jstat > 0) {
		cerr << "WARNING: warning occured inside call to SISL routine s1011. \n" 
		     << endl;
	    }

	    // write the result to file
	    writeGoCurve(result_curve, os);

	    freeCurve(result_curve);

	}

	os.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
