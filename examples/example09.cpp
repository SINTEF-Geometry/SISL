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
    string OUT_FILE_SURFACE = "example9_surf.g2";
    string OUT_FILE_POINTS = "example9_points.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program generates a series of surfaces that \n"
    "each interpolate an array of (internally defined) points in \n"
    "3D space.  The results will be written to file in Go-format.\n"
    "The generated surfaces will be written to the file \n"
    + OUT_FILE_SURFACE + "'.  The interpolated points will be \n"
    "written to the file '" + OUT_FILE_POINTS + "'.  The \n"
    "interpolating routine used is s1534.  When inspecting the \n"
    "result, note the effect of the varying order of the surfaces. \n"
    "Note also the 'artificial' bumps generated due to the \n"
    "interpolation criterion.\n\n";

    const int num_points_u = 5; // number of points in the u parameter direction
    const int num_points_v = 5; // number of points in the v parameter direction
    
    double points[] = 
	{
	 0, 0, 0,      1, 0, 0,      2, 0, 0,      3, 0, 0,      4, 0, 0,
	 0, 1, 0,      1, 1, 0,      2, 1, 0,      3, 1, 0,      4, 1, 0,
	 0, 2, 0,      1, 2, 0,      2, 2, 1,      3, 2, 0,      4, 2, 0,
	 0, 3, 0,      1, 3, 0,      2, 3, 0,      3, 3, 0,      4, 3, 0,
	 0, 4, 0,      1, 4, 0,      2, 4, 0,      3, 4, 0,      4, 4, 0
	};

    double u_par[] = {0, 1, 2, 3, 4}; // point parametrization in u-direction
    double v_par[] = {0, 1, 2, 3, 4}; // point parametrization in v-direction

    const int dim = 3; // dimension of the space we are working in
    
    const int num_surf = 4;

    const int order_u[] = {2, 3, 4, 4};
    const int order_v[] = {2, 3, 4, 2};

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
	ofstream os_pts(OUT_FILE_POINTS.c_str());
	if (!os_surf || !os_pts) {
	    throw runtime_error("Unable to open output file.");
	}

	// generating interpolating surface
	for (int i = 0; i < num_surf; ++i) {

	    SISLSurf* result_surf = 0;
	    int jstat = 0;

	    s1537(points,       // pointer to the array of points to interpolate
		  num_points_u, // number of interpolating points along the 'u' parameter
		  num_points_v, // number of interpolating points along the 'v' parameter
		  dim,          // dimension of the Euclidean space
		  u_par,        // pointer to the 'u' parameter values of the points
		  v_par,        // pointer to the 'v' parameter values of the points
		  0,            // no additional condition along edge 1
		  0,            // no additional condition along edge 2
		  0,            // no additional condition along edge 3
		  0,            // no additional condition along edge 4
		  order_u[i],   // the order of the generated surface in the 'u' parameter
		  order_v[i],   // the order of the generated surface in the 'v' parameter
		  1,            // open surface in the u direction
		  1,            // open surface in the v direction 
		  &result_surf, // the generated surface
		  &jstat);      // status variable
	    if (jstat < 0) {
		throw runtime_error("Error occured inside call to SISL routine s1537.");
	    } else if (jstat > 0) {
		cerr << "WARNING: warning occured inside call to SISL routine s1537. \n" 
		     << endl;
	    }
	    
	    // writing surface and interpolated points to file
	    writeGoSurface(result_surf, os_surf);
	    
	    freeSurf(result_surf);
	}

	writeGoPoints(num_points_u * num_points_v, points, os_pts);

	// cleaning up
	os_surf.close();
	os_pts.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
