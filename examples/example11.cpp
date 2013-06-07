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
    string OUT_FILE_SURFACE = "example11_surf.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program generates a cylindrical surface and saves it \n"
    "to the file '" + OUT_FILE_SURFACE + "'.\n\n";
    //==========================================================

    double bottom_pos[] = {10, 1, -1};
    double bottom_axis[] = {5, 0, 0};
    double ellipse_ratio = 1.2;
    double axis_dir[] = {0, 0, 1};
    double height = 35;

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {
	ofstream os_surf(OUT_FILE_SURFACE.c_str());
	if (!os_surf) {
	    throw runtime_error("Unable to open output file.");
	}

	SISLSurf* result_surf = 0;
	int jstat = 0;
	
	// generating the surface
	s1021(bottom_pos,      // bottom position of the cylinder
	      bottom_axis,     // the axis of the base of the cylinder
	      ellipse_ratio,   // the ellipticalness of the base of the cylinder
	      axis_dir,        // direction of the main axis
	      height,          // height of the cylinder
	      &result_surf,    // resulting surface
	      &jstat);         // status variable
	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1024.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1024. \n" 
		 << endl;
	}

	// writing surface to file
	writeGoSurface(result_surf, os_surf);

	// cleaning up
	freeSurf(result_surf);
	os_surf.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};

