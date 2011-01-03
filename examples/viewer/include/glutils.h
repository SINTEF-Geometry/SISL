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

#ifndef GLUTILS_H_INCLUDED


#ifndef _MSC_VER // 041103
#  include <GL/glx.h>
#endif
#include "GL/glut.h"

#include <vector>
using std::vector;

#include "jonvec.h"


//
// 991117: Note!!! glGetError shouldn't be called between glBegin/End!!!
//

//
// If gl-error, write out a message and abort execution.
// If no error, this should be fast.
// If you want to fool yourself, use the macro with function syntax!
//
void assert_gl_m(int line, char *file, const bool do_exit=true);
void assert_gl_dummy_and_empty(void);


//
// No semicolon; we won't allow anybody to destroy the source formatting...
//
#define ASSERT_GL assert_gl_m(__LINE__, __FILE__)
#define ASSERT_GL_DONT_EXIT assert_gl_m(__LINE__, __FILE__, false)

//
// This must be called with parantheses, function-like...
//
#define assert_gl assert_gl_m(__LINE__, __FILE__); assert_gl_dummy_and_empty

#ifndef _MSC_VER // 041103
void list_FBConfigs(GLXFBConfig *config, int nelements);
#endif

void draw_sphere(const vector3t<float> &pos, const float r,
		 const int n, const int slot,
		 const vector3t<float> &col,
		 const bool wiremode);






#define GLUTILS_H_INCLUDED
#endif
