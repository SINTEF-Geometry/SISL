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

#ifndef _GOREADWRITE_H
#define _GOREADWRITE_H

#include <iostream>
#include <vector>

// forward declarations
struct SISLCurve;
struct SISLSurf;

//===========================================================================
// For the two following read-functions, the memory for the SISL objects will be
// dynamically allocated inside the function, and a pointer is passed to the user.
// It is the user's responsibility to free the memory afterwards.  
// IMPORTANT:  Since the objects are allocated within SISL (which is written in 
// pure C), the C memory allocation functions are used.  To 'free' a SISLCurve 
// pointed to by 'c', you should call 'freeCurve(c)'.  The command to free
// SISLSurf objects is 'freeSurf()'.
//===========================================================================

// Read a curve in Go-format from stream and generate a SISLCurve object from it.
SISLCurve* readGoCurve(std::istream& go_stream);

// Read a surface in Go-format from stream and generate a SISLSurf object from it.
SISLSurf* readGoSurface(std::istream& go_stream);


//===========================================================================
// The following four functions are used to write SISL objects to a stream in
// Go-format.
//===========================================================================

// Write a SISLCurve object to stream using the Go-format
void writeGoCurve(SISLCurve* curve, std::ostream& go_stream);

// Write the SISLSurf object to stream using the Go-format.
void writeGoSurface(SISLSurf* surf, std::ostream& go_stream);

//===========================================================================
// The following is just an access function to read/write a set of 3D points
// on the Go format.  It is unrelated to SISL, but useful when you want to
// visualise 3D-points with the Go-compatible viewer.
//===========================================================================

// The 3D points whose coords are pointed to by 'coords' will be written in 
// Go-format to the stream.  The number of points is indicated by 'num_points'.
void writeGoPoints(int num_points, double* coords, std::ostream& go_stream);

// Read 3D coordinates grom a point cloud in Go-format.  The points' coordinates 
// will be filled into the vector 'coords' (whose length, divided by 3, is the 
// total number of points).  The points are read from the 'go_stream' stream.
void readGoPoints(std::vector<double>& coords, std::istream& go_stream);

#endif // _GOREADWRITE_H

