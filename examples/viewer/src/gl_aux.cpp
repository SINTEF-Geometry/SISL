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

#include <math.h>
#include <stdlib.h>
#include <cstring>

#include "aux2.h"
#include "gl_aux.h"



#ifndef PI
#  define PI 3.1415926536	// Where should it *really* be defined?
#endif


const int predefined_colours=24;
const double predef_col_table[3*predefined_colours]
={1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0,
  1.0, 1.0, 0.0,
  0.0, 1.0, 1.0,
  1.0, 0.0, 1.0,
  0.7, 0.0, 0.0,
  0.0, 0.7, 0.0,
  0.0, 0.0, 0.7,
  0.7, 0.7, 0.0,
  0.0, 0.7, 0.7,
  0.7, 0.0, 0.7,
  0.3, 0.8, 0.8,
  0.8, 0.3, 0.8,
  0.8, 0.8, 0.3,
  0.3, 0.3, 0.8,
  0.8, 0.3, 0.3,
  0.3, 0.8, 0.3,
  0.5, 0.8, 0.8,
  0.8, 0.5, 0.8,
  0.8, 0.8, 0.5,
  0.5, 0.5, 0.8,
  0.8, 0.5, 0.5,
  0.5, 0.8, 0.5};





//----------------------------------------------------------------------
//
// Drawing a cylinder, or possibly a cone.
//
//----------------------------------------------------------------------

void draw_cylinder(double x0, double y0, double z0,
		   double x1, double y1, double z1,
		   double radius, double radius2, int n)
{
  int i;
  double y, z, r;
  
  glPushMatrix();
  glTranslatef(x0, y0, z0);
  /*
    Now, we have to move P=(x1-x0, y1-y0, z1-z0) to (r, 0, 0) where
    r=length(x1-x0, y1-y0, z1-z0). First, rotate P to Q in the xy-plane:
  */
  /* Hvorfor blir det ikke riktig med neg. rotasjon her? */
  glRotatef(atan2(x1-x0, z1-z0)/PI*180.0-90.0,
	    0.0, 1.0, 0.0);
/*  printf("%f %f %f\t",
	 x1-x0,
	 z1-z0, 
	 atan2(x1-x0, z1-z0)/PI*180.0);*/
  /*
    Now, we have to move Q=(sqrt((x1-x0)^2+(z1-z0)^2), y1-y0, 0) to
    (r, 0, 0).
  */
  r=sqrt((x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
  /* Hvorfor blir det ikke riktig med neg. rotasjon her? */
  glRotatef(atan2(y1-y0, sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)))/PI*180.0,
	    0.0, 0.0, 1.0);
/*  printf("%f %f %f\n",
	 y1-y0,
	 sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)),
	 atan2(y1-y0, sqrt((x1-x0)*(x1-x0)+(z1-z0)*(z1-z0)))/PI*180.0);*/
  /*
    The next step is to draw the cylinder. Or cone. Normals not perfect for
    the cone. Impossible at pointy end with tri-fan?
  */
  if (radius==radius2)
    {
      glBegin(GL_QUAD_STRIP);
      for (i=0; i<=n; i++)
	{
	  y=radius*cos((2.0*PI*i)/n);
	  z=radius*sin((2.0*PI*i)/n);
	  glNormal3f(0.0, y, z);
	  glVertex3f(0.0, y, z);
	  glNormal3f(0.0, y, z);
	  glVertex3f(  r, y, z);
	}
      glEnd();
    }
  else
    {
      glBegin(GL_TRIANGLE_FAN);
      glNormal3f(1.0, 0.0, 0.0);
      glVertex3f(r, 0.0, 0.0);
      for (i=0; i<=n; i++)
	{
	  y=radius*cos((2.0*PI*i)/n);
	  z=radius*sin((2.0*PI*i)/n);
	  glNormal3f(0.0, y, z);
	  glVertex3f(0.0, y, z);
	}
      glEnd();
    }
  /*
    Finally, we restore the transformation stack.
  */
  glPopMatrix();
}






//----------------------------------------------------------------------
//
// Drawing a set of axes.
//
// 021203: Actually, this is newer. Fix!!!!
//
//----------------------------------------------------------------------

void draw_gl_axes_old(int n, double r, double radius, double rim, double l)
{
  // int n=10;
  // double r=0.7;
  // double radius=0.01, rim=0.04, l=0.1;
  
  const double vertex_red=1.0, vertex_green=1.0, vertex_blue=1.0;
  
  glColor3f(vertex_red, vertex_green, vertex_blue);
  
  draw_cylinder(0.0, 0.0, -r, 0.0, 0.0, r, radius, radius, n);
  draw_cylinder(0.0, 0.0, r-0.1, 0.0, 0.0, r+0.3, rim, 0.0, n);
  draw_cylinder(-0.5*l, l, r+0.3, 0.5*l, l, r+0.3,
		radius, radius, n);
  draw_cylinder(-0.5*l, l*2.0, r+0.3, 0.5*l, l*2.0, r+0.3,
		radius, radius, n);
  draw_cylinder(-0.5*l, l, r+0.3, 0.5*l, l*2.0, r+0.3,
		radius, radius, n);
  
  draw_cylinder(-r, 0.0, 0.0, r, 0.0, 0.0, radius, radius, n);
  draw_cylinder(r-0.1, 0.0, 0.0, r+0.3, 0.0, 0.0, rim, 0.0, n);
  draw_cylinder(r+0.3, l, 0.5*l, r+0.3, 2*l, -0.5*l,
		radius, radius, n);
  draw_cylinder(r+0.3, l, -0.5*l, r+0.3, 2*l, 0.5*l,
		radius, radius, n);
  
  draw_cylinder(0.0, -r, 0.0, 0.0, r, 0.0, radius, radius, n);
  draw_cylinder(0.0, r-0.1, 0.0, 0.0, r+0.3, 0.0, rim, 0.0, n);
  draw_cylinder(0.0, r+0.3, 2.0*l, 0.0, r+0.3, 1.5*l,
		radius, radius, n);
  draw_cylinder(0.0, r+0.3, 1.5*l, -0.5*l, r+0.3, l,
		radius, radius, n);
  draw_cylinder(0.0, r+0.3, 1.5*l, 0.5*l, r+0.3, l,
		radius, radius, n);
  
  glColor3f(1.0, 1.0, 1.0);
}






//----------------------------------------------------------------------
//
// Draw gridlines to visualize cells.
//
//----------------------------------------------------------------------

void draw_grid(const int n1, const int n2, const int n3)
{
  int i, j, k;
  
  glBegin(GL_LINES);
  for (i=0; i<n1; i++)
    for (j=0; j<n2; j++)
      {
	const double x=2.0*(i/(n1-1.0)-0.5);
	const double y=2.0*(j/(n2-1.0)-0.5);

	glColor3f(  1.0, 0.0,  0.0);
	glVertex3f(   x,   y, -1.0);
	glColor3f(  1.0, 0.0,  0.0);
	glVertex3f(   x,   y,  1.0);
      }
  for (j=0; j<n2; j++)
    for (k=0; k<n3; k++)
      {
	const double y=2.0*(j/(n2-1.0)-0.5);
	const double z=2.0*(k/(n3-1.0)-0.5);

	glColor3f(  0.0, 1.0, 0.0);
	glVertex3f(-1.0,   y,   z);
	glColor3f(  0.0, 1.0, 0.0);
	glVertex3f( 1.0,   y,   z);
      }
  for (k=0; k<n3; k++)
    for (i=0; i<n1; i++)
      {
	const double z=2.0*(k/(n3-1.0)-0.5);
	const double x=2.0*(i/(n1-1.0)-0.5);

	glColor3f(  0.0,  0.0, 1.0);
	glVertex3f(   x, -1.0,   z);
	glColor3f(  0.0,  0.0, 1.0);
	glVertex3f(   x,  1.0,   z);
      }
  glEnd();
}






//
// 020505: Adding transparent planes.
//

void draw_grid_planes(const int n1, const int n2, const int n3)
{
  int i, j, k;
  
  glDisable(GL_LIGHTING); // must be turned on by the caller...

  glBlendFunc(GL_ONE, GL_ONE);
  glEnable(GL_BLEND);
  glDisable(GL_DEPTH_TEST);
  glBegin(GL_QUADS);
  for (k=0; k<n3; k++)
    {
      const double z=2.0*(k/(n3-1.0)-0.5);
      
      glColor3f(  0.0,  0.0, 0.2);
      glVertex3f(-1.0, -1.0,   z);
      glVertex3f( 1.0, -1.0,   z);
      glVertex3f( 1.0,  1.0,   z);
      glVertex3f(-1.0,  1.0,   z);
    }
  for (j=0; j<n2; j++)
    {
      const double y=2.0*(j/(n2-1.0)-0.5);
      
      glColor3f(  0.0,  0.2,  0.0);
      glVertex3f(-1.0,    y, -1.0);
      glVertex3f( 1.0,    y, -1.0);
      glVertex3f( 1.0,    y,  1.0);
      glVertex3f(-1.0,    y,  1.0);
    }
  for (i=0; i<n1; i++)
    {
      const double x=2.0*(i/(n1-1.0)-0.5);
      
      glColor3f(  0.2,  0.0,  0.0);
      glVertex3f(   x, -1.0, -1.0);
      glVertex3f(   x,  1.0, -1.0);
      glVertex3f(   x,  1.0,  1.0);
      glVertex3f(   x, -1.0,  1.0);
    }
  glEnd();
  glEnable(GL_DEPTH_TEST);
  glDisable(GL_BLEND);
}






//----------------------------------------------------------------------
//
// Initializing glut and gl.
//
//----------------------------------------------------------------------

void gl_init(int argc, char *argv[],
	     const int xs, const int ys,
	     const int x0, const int y0,
	     const int doubleBuffer,
	     const int two_sided,
	     const int lighting,
	     const int normalize,
	     const int smooth,
	     const double xtrans,
	     const double ytrans,
	     const double ztrans,
	     const double xscale,
	     const double yscale,
	     const double zscale,
	     int &tx,
	     int &ty,
	     const int texture_mode,
	     const char * const texfile /* =NULL */)
{
  int i;
  const int use_accum=0;
  
  //
  // 020430: NB! Note that the really confusing name 'GLUT_RGBA' hides the
  //         fact that this setting does not induce an alpha buffer being set
  //         up! One must in fact specify GLUT_ALPHA explicitly!
  //         (But I'm not sure how important this is... Perhaps this alpha is
  //         not actually needed for the *window* for OpenGL for instance to
  //         be able to do alpha blending?!)
  // 020501: Seems that we get accum buffer even without specifying it?!
  // 020502: GLUT_ALPHA seems not to be available on smaug.
  //         But we don't actually need it for alpha-blending, anyway. For
  //         that we only need alpha for the textures.
  //         Ok, it's there when we go to 24 bpp mode in X.
  //
#ifndef QTMODE
  glutInitDisplayMode(GLUT_RGBA | // GLUT_ALPHA | 
		      (use_accum ? GLUT_ACCUM : 0) |
		      ((doubleBuffer) ? GLUT_DOUBLE : GLUT_SINGLE) |
		      GLUT_DEPTH);
  glutInitWindowSize(xs, ys);
  glutInitWindowPosition(x0, y0);
  glutInit(&argc, argv);
  glutCreateWindow(argv[0]);
#endif

  glClearColor(0.0, 0.0, 0.0, 0.0);
  if (use_accum)
    glClearAccum(0.0, 0.0, 0.0, 0.0);
  if (smooth)
    glShadeModel(GL_SMOOTH);
  else
    glShadeModel(GL_FLAT);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_LESS);			/* GL_LESS is default. */
/*  glDepthRange(100.0, -100.0);*/	/* What is default value? */
/*  glClearDepth(0.0); */		/* What is default value? */

  /* Initialize materials */

  {
    static float ambient[] = {0.1, 0.1, 0.1, 1.0};
    static float diffuse[] = {0.4, 0.4, 0.4, 1.0};

    static float position0[] = { 5.0, 5.0, 200.0, 0.0};
    static float position1[] = {-5.0, -5.0, 200.0, 0.0};

    static float front_mat_shininess[] = {50.0};
    // static float front_mat_specular[] = {0.9, 0.9, 0.9, 1.0};
    static float front_mat_specular[] = {0.0, 0.0, 0.0, 1.0};
    static float front_mat_diffuse[] = {0.7, 0.7, 0.7, 1.0};
    // static float front_mat_diffuse[] = {0, 0, 0, 1.0};
    static float front_mat_ambient[] = {0.3, 0.3, 0.3, 1.0};
    static float front_mat_emission[] = {0.3, 0.0, 0.0, 1.0};

    static float back_mat_shininess[] = {128.0};
    // static float back_mat_specular[] = {0.4, 0.4, 0.4, 1.0};
    static float back_mat_specular[] = {0.0, 0.0, 0.0, 1.0};
    static float back_mat_diffuse[] = {0.7, 0.7, 0.7, 1.0};
    static float back_mat_ambient[] = {0.3, 0.3, 0.3, 1.0};
    static float back_mat_emission[] = {0.0, 0.3, 0.0, 1.0};
  
    static float lmodel_ambient[] = {1.0, 1.0, 1.0, 1.0};
    static float lmodel_twoside[] = {GL_TRUE};

    if (two_sided)
      lmodel_twoside[0]=GL_TRUE;
    else
      lmodel_twoside[0]=GL_FALSE;
    
    glLightfv(GL_LIGHT0, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT0, GL_POSITION, position0);
    
    glLightfv(GL_LIGHT1, GL_AMBIENT, ambient);
    glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
    glLightfv(GL_LIGHT1, GL_POSITION, position1);

    glLightModelfv(GL_LIGHT_MODEL_AMBIENT, lmodel_ambient);
    glLightModelfv(GL_LIGHT_MODEL_TWO_SIDE, lmodel_twoside);

    // glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);

    glMaterialfv(GL_FRONT, GL_SHININESS, front_mat_shininess);
    glMaterialfv(GL_FRONT, GL_SPECULAR,  front_mat_specular);
    glMaterialfv(GL_FRONT, GL_DIFFUSE,   front_mat_diffuse);
    glMaterialfv(GL_FRONT, GL_AMBIENT,   front_mat_ambient);
    glMaterialfv(GL_FRONT, GL_EMISSION,  front_mat_emission);

    glMaterialfv(GL_BACK, GL_SHININESS, back_mat_shininess);
    glMaterialfv(GL_BACK, GL_SPECULAR,  back_mat_specular);
    glMaterialfv(GL_BACK, GL_DIFFUSE,   back_mat_diffuse);
    glMaterialfv(GL_BACK, GL_AMBIENT,   back_mat_ambient);
    glMaterialfv(GL_BACK, GL_EMISSION,  back_mat_emission);

    glEnable(GL_LIGHT0);
    glEnable(GL_LIGHT1); 
    if (lighting)
      glEnable(GL_LIGHTING);
    else
      glDisable(GL_LIGHTING);
    if (normalize)
      glEnable(GL_NORMALIZE);
    else
      glDisable(GL_NORMALIZE);
  }
  
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
/*  glFrustum( -1.0, 1.0, -1.0, 1.0, 5, 25 ); */
/*  glFrustum( -1.0, 1.0, -1.0, 1.0, 2, 4 ); */
/*  glFrustum( -1.0, 1.0, -1.0, 1.0, 100, 102 ); */
/*  glFrustum( -1.0, 1.0, -1.0, 1.0, 5, 7 );*/
  /*
    Note that the far clipping plane distance doesn't influence on
    the projection, but that the distance between near and far should be
    as close to one as possible to best utilise the z-buffer. Or something
    like this...
  */
  //
  // 000304: Was near=5, near clipping plane to far from viewer.
  //         6-0.5*sqrt(3) should make a 1x1x1 cube never intersect the near
  //         plane regardless of rotation... Making it even closer will
  //         make it possible to expand things more before being clipped...
  //
  //glFrustum(-1.0, 1.0, -1.0, 1.0, 6-0.5*sqrt(3), 20);
  
  //
  // 020214: This one seems to be rather good, we can zoom in quite a lot
  //         without stuff hitting the near clipping plane.
  //
  glFrustum(-1.0, 1.0, -1.0, 1.0, 3, 20); // cmt out 020201

  // 020410 experimenting
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 3, 20);

  //
  // 020214: The next one let's us see inside more easily, because the near
  //         clipping plane is further away from us.
  //
  //glFrustum(-1.0, 1.0, -1.0, 1.0, 5, 20);

  //
  // 021206: Trying to get the near clipping plane closer to us.
  //         I've lost track of the actual numbers here...
  //         Should be dealt with...
  //
  glLoadIdentity();
  glFrustum(-1.0, 1.0, -1.0, 1.0, 2.0, 20);

  //
  // 020627: Trying to reduce the "perspective" effect.
  //
  //glFrustum(-1.0, 1.0, -1.0, 1.0, 5+100, 20+100);
  
  glMatrixMode(GL_MODELVIEW);
  glLoadIdentity();
/*  glTranslated( 0.0, 0.0, -6.0 ); */
/*  glTranslated( 0.0, 0.0, -3.0 ); */
/*  glTranslated( 0.0, 0.0, -101.0 ); */
  glTranslated(xtrans, ytrans, ztrans);
  glScaled(xscale, yscale, zscale);

  /*
    Why did I do this? To clear both buffers? Is this necessary here?
  */
  for (i=0; i<2; i++)
    {
      glDrawBuffer(GL_BACK);
      glReadBuffer(GL_BACK);
      glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
      glClear(GL_COLOR_BUFFER_BIT | (use_accum ? GL_ACCUM_BUFFER_BIT : 0));
#ifndef QTMODE
      glutSwapBuffers();
#endif
    }

  //
  // Note that after the texture has been inserted into the GL engine, it
  // is lost/deallocated from the application.
  // (But (for some reason) the size is returned...)
  //
  unsigned char *texture=NULL;
  if (texfile==NULL)
    {
      if (tx!=ty)
	CRIT_ERR(printf("tx=%d, ty=%d, not implemented.\n", tx, ty));
      
      if ((texture=new unsigned char[3*SQR(tx)])==NULL)
	CRIT_ERR(puts("Couldn't allocate space for texture."));
      
      int i;
      const unsigned char a=255, b=0;
      
      for (i=0; i<3*SQR(tx); i++)
	texture[i]=a;
      for (i=0; i<tx; i++)
	{
	  texture[3*(i*tx+0)+0]=b;
	  texture[3*(i*tx+0)+1]=b;
	  texture[3*(i*tx+0)+2]=b;
	  texture[3*(i*tx+(tx-1))+0]=b;
	  texture[3*(i*tx+(tx-1))+1]=b;
	  texture[3*(i*tx+(tx-1))+2]=b;
	  texture[3*(0*tx+i)+0]=b;
	  texture[3*(0*tx+i)+1]=b;
	  texture[3*(0*tx+i)+2]=b;
	  texture[3*((tx-1)*tx+i)+0]=b;
	  texture[3*((tx-1)*tx+i)+1]=b;
	  texture[3*((tx-1)*tx+i)+2]=b;
	}
    }
  else
    //
    // Read texture from file...
    //
    {
      puts("texture code not compiled in!");
      exit(0);
#if 0
      unsigned char *texture_tmp=read_ppm_file(texfile, &tx, &ty);
      puts("Leste pgm-fil.");
      if (tx!=ty)
	CRIT_ERR(puts("Not implemented."));
      if ((tx!=1) && (tx!=2) && (tx!=4) && (tx!=8) && (tx!=16) && (tx!=32) && 
	  (tx!=64) && (tx!=128) && (tx!=256) && (tx!=512) && (tx!=1024) &&
	  (tx!=2048) && (tx!=4096))
	CRIT_ERR(puts("Texture size not ok."));

      printf("Texture size is %dx%d.\n", tx, ty);
      
      if ((texture=new unsigned char[3*SQR(tx)])==NULL)
	CRIT_ERR(puts("Couldn't allocate space for texture."));

      int i;
     
      for (i=0; i<SQR(tx); i++)
	{
	  texture[3*i+0]=texture_tmp[i];
	  texture[3*i+1]=texture_tmp[i];
	  texture[3*i+2]=texture_tmp[i];
	  // printf("%5d - %3d %3d %3d\n", i,
	  // texture[i+0], texture[i+1], texture[i+2]);
	}

      delete texture_tmp;
#endif
    }
  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tx, ty,
	       0, GL_RGB, GL_UNSIGNED_BYTE, texture);
  glTexEnvi(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, texture_mode);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  
  delete texture;

  //glPolygonOffset(0.0, 20.0);
  //
  // 020430: PolygonOffset(factor, units) produces an offset o, where
  //         o = max depth slope * factor + impl. dep. const. r * units.
  //         0.0, 20.0 worked just perfect for GeForce 2 and  Matrox Mystique,
  //         but not for GeForce 3.
  //         Perhaps it was strange that it worked? Or maybe it enough now
  //         just use a factor slightly greater than 0.0?
  //         A good idea would perhaps be to make a tool to adjust these
  //         empirically during runtime?!
  //         Seems that it is not enough to just increase factor slightly,
  //         perhaps the GF3 has another definition for 'units'? Increasing
  //         it somewhat seems to make the problem go away...
  //         Don't want to increase it too much though.
  //         Perhaps a way to do this: Increase the 'factor' until very
  //         slanted surfaces turns out ok, and increase 'units' until
  //         flat (parallel to near plane) surfaces are ok?! Will do it
  //         like this but opposite order.
  //         Hmm!!! (1.0, 1.0) seems to be just fine!
  // 
  glPolygonOffset(1.0, 1.0);
  glEnable(GL_POLYGON_OFFSET_POINT);
  glEnable(GL_POLYGON_OFFSET_LINE);
  glEnable(GL_POLYGON_OFFSET_FILL);

}






void transpose_matrix(double * const d)
{
  int i, j;
  
  for (i=0; i<4; i++)
    for (j=0; j<i; j++)
      {
	double tmp=d[i*4+j];
	d[i*4+j]=d[j*4+i];
	d[j*4+i]=tmp;
      }
}






void print_gl_matrix(const int m)
{
  int i;
  double p[16], p2[16];
  
  switch (m)
    {
    case GL_PROJECTION_MATRIX:

    case GL_MODELVIEW_MATRIX:
      glGetDoublev((GLenum)m, p);
      break;

#if 0
    case GL_MODELVIEW_MATRIX + GL_PROJECTION_MATRIX:
      glMatrixMode(GL_MODELVIEW);
      glPushMatrix();
      glGetDoublev(GL_PROJECTION_MATRIX, p);
      glMultMatrixd(p);
      glGetDoublev(GL_MODELVIEW_MATRIX, p);
      glPopMatrix();
      break;
#endif

    case GL_MODELVIEW_MATRIX + GL_PROJECTION_MATRIX:
      glMatrixMode(GL_MODELVIEW);
      glGetDoublev(GL_MODELVIEW_MATRIX, p2);
      glPushMatrix();
      glLoadIdentity();
      glGetDoublev(GL_PROJECTION_MATRIX, p);
      glMultMatrixd(p);
      glMultMatrixd(p2);
      glGetDoublev(GL_MODELVIEW_MATRIX, p);
      glPopMatrix();
      break;
    }

  //
  // 000320: BUG was here! GL uses column-order...
  //
  transpose_matrix(p);
  
  printf("{");
  for (i=0; i<15; )
    {
      printf("%g, ", p[i]);
      i++;
      if (i%4==0)
	printf("\n");
    }
  printf("%g}\n\n", p[15]);

  printf("For Mathematica:\n\n");
  printf("{");
  for (i=0; i<15; )
    {
      if (i%4==0)
	printf("{");
      printf("%f", p[i]);
      i++;
      if (i%4==0)
	printf("}, ");
      else
	printf(", ");
    }
  printf("%f}}\n\n", p[15]);
 

}







//----------------------------------------------------------------------
//
// This callback is called when the window is reshaped.
//
//----------------------------------------------------------------------

void reshape_window(int width, int height)
{
  glViewport(0, 0, (GLint)width, (GLint)height);
/*
  glutPostRedisplay();
*/
}






//----------------------------------------------------------------------
//
// Writing and reading the GL transformation matrices to and from
// files.
//
//----------------------------------------------------------------------

void write_gl_matrices(FILE *f)
{
  int i, current_mode;
  double p[16];
  
  glGetIntegerv(GL_MATRIX_MODE, &current_mode);

  fprintf(f, "GL_PROJECTION_MATRIX\n");
  glGetDoublev(GL_PROJECTION_MATRIX, p);
  for (i=0; i<16; i++)
    fprintf(f, "%.15e\n", p[i]);

  fprintf(f, "GL_MODELVIEW_MATRIX\n");
  glGetDoublev(GL_MODELVIEW_MATRIX, p);
  for (i=0; i<16; i++)
    fprintf(f, "%.15e\n", p[i]);

  glMatrixMode(current_mode);
}

void read_gl_matrices(FILE *f)
{
  int i, current_mode;
  double p[16];
  char s[1000];
  
  glGetIntegerv(GL_MATRIX_MODE, &current_mode);

  fgets(s, 1000, f);
  if (strcmp(s, "GL_PROJECTION_MATRIX\n")!=0)
    CRIT_ERR(printf("Corrupt file contents? Read '%s'.\n", s));
  for (i=0; i<16; i++)
    fscanf(f, "%lf\n", p+i);
  glMatrixMode(GL_PROJECTION_MATRIX);
  glLoadMatrixd(p);

  fgets(s, 1000, f);
  if (strcmp(s, "GL_MODELVIEW_MATRIX\n")!=0)
    CRIT_ERR(printf("Corrupt file contents? Read '%s'.\n", s));
  for (i=0; i<16; i++)
    fscanf(f, "%lf\n", p+i);
  glMatrixMode(GL_MODELVIEW_MATRIX);
  glLoadMatrixd(p);

  glMatrixMode(current_mode);
}
