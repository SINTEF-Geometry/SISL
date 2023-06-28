/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include <iostream>
#include <fstream>
#include <string>
#include <stdexcept>
#include <cstdlib>

#include "sisl.h"
#include "GoReadWrite.h"

using namespace std;



//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
  ifstream is_sf1(vararg[1]);
  ifstream is_sf2(vararg[2]);
  double eps = atof(vararg[3]);
  ofstream os(vararg[4]);

  try {


    // reading surfaces
    SISLSurf* surf_1 = readGoSurface(is_sf1);
    SISLSurf* surf_2 = readGoSurface(is_sf2);

    // detecting (but not tracing out) intersection curves
    double epsco = 1.0e-15; // computational epsilon
    double epsge = 1.0e-5; // geometric tolerance
    int num_int_points = 0; // number of detected intersection points
    double* intpar_surf_1  = 0; // parameter values for the surface in the intersections
    double* intpar_surf_2 = 0; // parameter values for the curve in the intersections
    int num_int_curves = 0;   // number of intersection curves
    SISLIntcurve** intcurve = 0; // pointer to array of detected intersection curves
    int jstat = 0; // status variable

    // calculating topology of intersections
    s1859(surf_1,          // the first surface
	  surf_2,          // the second surface
	  epsco,           // computational resolution
	  eps,           // geometry resolution
	  &num_int_points, // number of single intersection points
	  &intpar_surf_1,  // pointer to array of parameter values for surface 1
	  &intpar_surf_2,  //               -"-                    for surface 2
	  &num_int_curves, // number of detected intersection curves
	  &intcurve,       // pointer to array of detected intersection curves.
	  &jstat);         // status variable
	
    if (jstat < 0) {
      throw runtime_error("Error occured inside call to SISL routine s1859.");
    } else if (jstat > 0) {
      cerr << "WARNING: warning occured inside call to SISL routine s1859. \n" 
	   << endl;
    }

    // In this example, we expect to detect two intersection curves, but no isolated
    // intersection points.

    cout << "Number of intersection points detected: " << num_int_points << endl;
    for (int i=0; i<num_int_points; ++i)
      cout << intpar_surf_1[2*i] << " " << intpar_surf_1[2*i+1] << " " << intpar_surf_2[2*i] << " " << intpar_surf_2[2*i+1] << endl;
    cout << "Number of intersection curves detected: " << num_int_curves << endl;
    for (int i=0; i<num_int_curves; ++i)
      {
	int npt = intcurve[i]->ipoint;
	cout << intcurve[i]->epar1[0] << " " << intcurve[i]->epar1[1] << " " ;
	cout << intcurve[i]->epar2[0] << " " << intcurve[i]->epar2[1] << " " ;
	// evaluating (tracing out) intersection curves and writing them to file
	cout << intcurve[i]->epar1[2*(npt-1)] << " " << intcurve[i]->epar1[2*npt-1] << " " ;
	cout << intcurve[i]->epar2[2*(npt-1)] << " " << intcurve[i]->epar2[2*npt-1] << endl;
      }
    for (int i = 0; i < num_int_curves; ++i) {
      s1310(surf_1,          // the first surface
	    surf_2,          // the second surface
	    intcurve[i],     // pointer to the intersection curve object 
	    epsge,           // geometric tolerance
	    double(0),       // maximum step size (ignored if <= 0)
	    1,               // make only 3D curve (no 2D curve in parametric domain)
	    0,               // don't draw the curve
	    &jstat);
      if (jstat < 0) {
	throw runtime_error("Error occured inside call to SISL routine s1310.");
      } else if (jstat == 3) {
	throw runtime_error("Iteration stopped due to singular point or degenerate "
			    "surface.");
      }
      writeGoCurve(intcurve[i]->pgeom, os);
    }

    // cleaning up
    if (surf_1) freeSurf(surf_1);
    if (surf_2) freeSurf(surf_2);
    is_sf1.close();
    is_sf2.close();
    os.close();
    if (intpar_surf_1) free(intpar_surf_1);
    if (intpar_surf_2) free(intpar_surf_2);
    if (num_int_curves > 0)
      freeIntcrvlist(intcurve, num_int_curves);

  } catch (exception& e) {
    cerr << "Exception thrown: " << e.what() << endl;
    return 0;
  }

  return 1;
};
