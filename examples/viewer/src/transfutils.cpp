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
