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
    string IN_FILE_CURVE_1 = "example1_curve.g2";
    string IN_FILE_CURVE_2 = "example2_curve.g2";
    string OUT_FILE_CURVE  = "example3_curve.g2";

    string DESCRIPTION = 
    "This program will create a 'blend curve' to connect the \n"
    "endpoints of the two curves generated in the two previous \n"
    "example programs.  The routine used is s1606. \n"
    "Input: " + IN_FILE_CURVE_1 + " and " + IN_FILE_CURVE_2 + "\n"
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
	ifstream stream_1(IN_FILE_CURVE_1.c_str());
	ifstream stream_2(IN_FILE_CURVE_2.c_str());
	if (!stream_1 || !stream_2) {
	    string error_message = 
		"Unable to open input files: " + IN_FILE_CURVE_1 +
		" and " + IN_FILE_CURVE_2 + ".  Are you sure you have run the "
		"two previous sample programs?";
	    throw runtime_error(error_message.c_str());
	}

	SISLCurve* c1 = readGoCurve(stream_1);
	SISLCurve* c2 = readGoCurve(stream_2);

	double epsge = 1.0e-5; // geometric precision
	int blendtype = 0; // generate polynomial segment
	int dim = 3;
	int order = 4;
	double c1_endpoint[3]; // endpoint of curve 1 (must be calculated)
	double c2_endpoint[3]; // endpoint of curve 2 (must be calculated)
	
	// The blend curve will extend from the endpoint of curve 1 to the
	// endpoint of curve 2.  As input, the SISL routine needs (approximate)
	// coordinates for these points on the curve.  We therefore preliminarly
	// need to evaluate these.  For this purpose, we use SISL routine s1227.
	
	// The end parameters of the curves' parametric domains can be found by 
	// looking at their knotvectors (pointed to by data member 'et').  
	// If the number of control points is 'n', then the knot numbered 'n'
	// would represent the end parameter (when counting from 0).  The number
	// of control points is indicated by the data member 'in'.
	double c1_endpar = c1->et[c1->in]; // end parameter of curve 1
	double c2_endpar = c2->et[c2->in]; // end parameter of curve 2
	int temp, jstat1, jstat2;
	
	// evaluating endpoint positions of both curves
	s1227(c1,          // input curve
	      0,           // evaluate position only (no derivatives)
	      c1_endpar,   // end parameter
	      &temp,       // indicates param. interval (not interesting for our purposes)
	      c1_endpoint, // this is what we want to calculate (3D position)
	      &jstat1);     // status variable (0 if everything all right)
	
	s1227(c2,          // input curve
	      0,           // evaluate position only (no derivatives)
	      c2_endpar,   // end parameter
	      &temp,       // indicates param. interval (not interesting for our purposes)
	      c2_endpoint, // this is what we want to calculate (3D position)
	      &jstat2);     // status variable (0 if everything all right)
	
	if (jstat1 < 0 || jstat2 < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1227.");
	} else if (jstat1 > 0 || jstat2 > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1227.";
	}
    
	// calculating blend curve
	SISLCurve* blend_curve = 0;
	int jstat;

	s1606(c1,           // the first input curve
	      c2,           // the second input curve
	      epsge,        // geometric tolerance
	      c1_endpoint,  // endpoint of curve 1 (geometric)
	      c2_endpoint,  // endpoint of curve 2 (geometric)
	      blendtype,    // type of blend curve (circle, conic, polynomial)
	      dim,          // dimension (3D)
	      order,        // order of generated spline curve
	      &blend_curve, // the generated curve
	      &jstat);      // status message
	     
	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1606.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1606.\n" << endl;
	}

	ofstream os(OUT_FILE_CURVE.c_str());
	if (!os) {
	    throw runtime_error("Unable to open output file.");
	}

	// write result to file
	writeGoCurve(blend_curve, os);
	
	// cleaning up
	freeCurve(blend_curve);
	freeCurve(c1);
	freeCurve(c2);
	os.close();
	stream_1.close();
	stream_2.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
