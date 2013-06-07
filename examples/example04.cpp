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
