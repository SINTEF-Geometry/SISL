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
#include <algorithm>
#include <stdexcept>

#include "sisl.h"
#include "GoReadWrite.h"

using namespace std;

namespace {
    string OUT_FILE_POINTS = "example14_points.g2";
    string OUT_FILE_CURVE = "example14_curve.g2";

    string DESCRIPTION = 
    //==========================================================
    "This program demonstrates one of the data reduction \n"
    "capabilities in SISL.  It generates a dense point set by \n"
    "sampling from a predefined curve, and then try to generate \n"
    "a curve that fits closely to to these samples, using as few \n"
    "control points as possible.  (In our case, since we know that \n"
    "the points were sampled from a simple curve, we already \n"
    "have knowledge about the ideal solution to this problem). \n"
    "The sampled points will be saved to the file '" + OUT_FILE_POINTS + "'\n"
    "and the generated curve will be saved to the file \n'"
    + OUT_FILE_CURVE + "'.\n\n";
    //==========================================================


    // information for the sample curve (a simple parabola)
    const int number = 3;
    const int order = 3;
    
    double coef[] = {0, 0, 0,
		     1, 1, 0,
		     2, 0, 0};
    
    double knots[] = {0, 0, 0, 1, 1, 1};
    
    const int num_samples = 500;

}; // end anonymous namespace 

//===========================================================================
int main(int avnum, char** vararg)
//===========================================================================
{
    cout << '\n' << vararg[0] << ":\n" << DESCRIPTION << endl;
    cout << "To proceed, press enter, or ^C to quit." << endl;
    getchar();

    try {
	ofstream os_pts(OUT_FILE_POINTS.c_str());
	ofstream os_cv(OUT_FILE_CURVE.c_str());
	if (!os_pts || !os_cv) {
	    throw runtime_error("Unable to open output file.");
	}
	
	// generating sample curve
	SISLCurve* sample_curve = 
	    newCurve(number, // number of control points
		     order,  // order of spline curve (degree + 1)
		     knots,  // pointer to knot vector (parametrization)
		     coef,   // pointer to coefficient vector (control points)
		     1,      // kind = polynomial B-spline curve
		     3,      // dimension
		     0);     // no copying of information, 'borrow' arrays
	if (!sample_curve) {
	    throw runtime_error("Error occured while generating curve.");
	}

	// sampling curve densely
	vector<double> sampled_point_coords(3 * num_samples);
	for (int i = 0; i < num_samples; ++i) {
	    double parameter = double(i)/(num_samples - 1);
	    int temp, jstat;

	    s1227(sample_curve, // the curve to sample from
		  0,            // calculate no derivatives
		  parameter,    // sample in this parameter
		  &temp,        // not used for our purposes (returns parameter interval)
		  &sampled_point_coords[3 * i], // sampled point written here
		  &jstat); // status
	    if (jstat < 0) {
		throw runtime_error("Error occured inside call to SISL routine s1227.");
	    } else if (jstat > 0) {
		cerr << "WARNING: warning occured inside call to SISL routine 1227."
		     << endl;
	    }
	}

	//generating approximating curve
	SISLCurve* result_curve = 0;
	const double tolerance = 1.0e-5; // pointwise tolerance
	const double afctol = 0.1; // number between 0 and 1 describing how tolerance
	                           // is to be distributed in two different steps of the
	                           // algorithm.
	const int max_iterations = 150;
	vector<double> tolerance_vec(num_samples, tolerance); // same tolerance for all samples
	vector<double> emxerr(num_samples);
	int jstat = 0;

	s1961(&sampled_point_coords[0], // array of points to be approximated
	      num_samples,              // number of sample points
	      3,                        // dimension of Euclidean space
	      2,                        // flag indicating uniform parametrization
	      NULL,                     // we do not use this argument
	      &tolerance_vec[0],        // the pointwise tolerance
	      0,                        // number of fixed derivatives at left
	      0,                        // number of fixed derivatives at right
	      1,                        // flag indicating that we want an open curve
	      afctol,                   // distribution of tolerance
	      max_iterations,           // maximum number of iterations in the routine
	      3,                        // polynomial order of approximation
	      &result_curve,            // pointer to the generated curve
	      &emxerr[0],               // reporting the max deviation in each point
	      &jstat);                  // status variable

	if (jstat < 0) {
	    throw runtime_error("Error occured inside call to SISL routine s1961.");
	} else if (jstat > 0) {
	    cerr << "WARNING: warning occured inside call to SISL routine s1961. \n" 
		 << endl;
	}
	
	cout << "Number of initial data points: " << num_samples << endl;
	cout << "Number of control points in final approximation: " << result_curve->in << endl;
	cout << "Max pointwise error is : " << *(max_element(emxerr.begin(), emxerr.end()));
	cout << endl << endl;

	// saving to file
	writeGoPoints(num_samples, &sampled_point_coords[0], os_pts);
	writeGoCurve(result_curve, os_cv);

	// cleaning up
	freeCurve(result_curve);
	os_pts.close();
	os_cv.close();

    } catch (exception& e) {
	cerr << "Exception thrown: " << e.what() << endl;
	return 0;
    }

    return 1;
};
