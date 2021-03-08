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

