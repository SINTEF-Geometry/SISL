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

#include <GL/glut.h>
#include <stdio.h>

#include "transfutils.h"






//
// 990831: These are used for storing "global" rotations etc. It's convenient
//         to do this with global variables, since the mouse callbacks can only
//         take parameters fixed by glut...
//
// 980913: ztrans -6 makes the near clipping plane "visible", because it's
//         at -5. ztrans -5 would fix this.
//
double xtrans=0.0, ytrans=0.0, ztrans=-6;
double xrot_eps=5.0, yrot_eps=5.0, zrot_eps=5.0;
double xscale=1.0, yscale=1.0, zscale=1.0;
double xrot, yrot, zrot; // These are initialized in 'Mouse'.






void rotate(double y_ang, double x_ang, double z_ang)
{
  float oldModelView[16];

  if (y_ang)
    {
      glGetFloatv(GL_MODELVIEW_MATRIX, oldModelView);
      glLoadIdentity();
      glTranslated(xtrans, ytrans, ztrans);
      glRotated(y_ang*yrot_eps, 0.0, 1.0, 0.0);
      glTranslated(-xtrans, -ytrans, -ztrans);
      glMultMatrixf(oldModelView);
    }
  
  if (x_ang)
    {
      glGetFloatv(GL_MODELVIEW_MATRIX, oldModelView);
      glLoadIdentity();
      glTranslated(xtrans, ytrans, ztrans);
      glRotated(x_ang*xrot_eps, 1.0, 0.0, 0.0);
      glTranslated(-xtrans, -ytrans, -ztrans);
      glMultMatrixf(oldModelView);
    }

  if (z_ang)
    {
      glGetFloatv(GL_MODELVIEW_MATRIX, oldModelView);
      glLoadIdentity();
      glTranslated(xtrans, ytrans, ztrans);
      glRotated(z_ang*zrot_eps, 0.0, 0.0, 1.0);
      glTranslated(-xtrans, -ytrans, -ztrans);
      glMultMatrixf(oldModelView);
    }
}






void translate(double x, double y, double z)
{
  float oldModelView[16];

  glGetFloatv(GL_MODELVIEW_MATRIX, oldModelView);
  glLoadIdentity();
  glTranslated(x, y, z);
  glTranslated(-xtrans, -ytrans, -ztrans);
  xtrans=x;
  ytrans=y;
  ztrans=z;
  glMultMatrixf(oldModelView);
}






void scale(double x, double y, double z)
{
  glScaled(x, y, z);
  xscale*=x;
  yscale*=y;
  zscale*=z;
}
