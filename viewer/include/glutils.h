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
