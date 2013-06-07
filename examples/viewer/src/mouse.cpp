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

#include <stdio.h>
#define GLUT_DISABLE_ATEXIT_HACK
#include <GL/glut.h>
#include <math.h>
#ifndef _MSC_VER
#  include <unistd.h>
#endif

#ifndef PI
#  define PI 3.1415926536
#endif
#include "mouse.h"
#include "transfutils.h"
#include "aux2.h"

//
// All casting to double in calls to sqrt are just because of aCC on tor.
// Doing this by using temporary variables instead...
//






//
// 990831: These are used for storage of flags, states, old mouse-
//         positions etc.
//
GLboolean allow_zrot=GL_TRUE;
int last_x=0, last_y=0, last_x0=0, last_y0=0;
bool mouse_movement=false; // Mouse-movement routines set this.






//
// 040326: An attempt to make a "draw cursor into framebuffer" subroutine,
//         in order to get the cursor into movie sequences produced with
//         'sfviewer1'.
//

float curx=0.0, cury=0.0, curz=0.0;

void draw_cursor(int x, int y)
{
  GLint viewport[4];
  GLdouble modelview[16];
  GLdouble projection[16];
  GLfloat winX, winY, winZ;
  GLdouble posX, posY, posZ;
  
  glGetDoublev( GL_MODELVIEW_MATRIX, modelview );
  glGetDoublev( GL_PROJECTION_MATRIX, projection );
  glGetIntegerv( GL_VIEWPORT, viewport );
  
  winX = (float)x;
  winY = (float)viewport[3] - (float)y;
  glReadPixels( x, int(winY), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ );
  
  gluUnProject(winX, winY, winZ,
	       modelview, projection, viewport, &posX, &posY, &posZ);

  curx=posX;
  cury=posY;
  curz=posZ;
}






int transversal_rotation(int x, int y, int last_xx, int last_yy)
{
  int viewp[4];

  glGetIntegerv(GL_VIEWPORT, viewp);
  if ((fabs(((x-last_xx)*(x-(viewp[0]+viewp[2])*0.5)+
	     (y-last_yy)*(y-(viewp[1]+viewp[3])*0.5))/
	    sqrt((double)(SQR(x-last_xx)+SQR(y-last_yy)))/
	    sqrt((double)(SQR(x-(viewp[0]+viewp[2])*0.5)+
			  SQR(y-(viewp[1]+viewp[3])*0.5)))) >
       cos(50.0/180.0*PI)) ||
      (sqrt((x-(viewp[0]+viewp[2])*0.5)*(x-(viewp[0]+viewp[2])*0.5)+
	    (y-(viewp[1]+viewp[3])*0.5)*(y-(viewp[1]+viewp[3])*0.5)) <
       0.2*sqrt((double)SQR(viewp[2]-viewp[0])+SQR(viewp[3]-viewp[1]))) ||
      (!allow_zrot))
    {
      /* puts("radielt");*/
      return 0;
    }
  else
    {
      /* puts("transversalt");*/
      return 1;
    }
  // draw_cursor(x, y);
}






//
// 990831: Would have liked to send other parameters too, but 
//         glut won't allow that...
//
void MouseRotate(int x, int y)
{
  int viewp[4];
  //const double sensitiveness=0.1;
  // 030129: Reducing sensitivity slightly.
  //const double sensitiveness=0.05;
  // 030302: Reducing sensitivity slightly more, for demo.
  const double sensitiveness=0.02;

  mouse_movement = mouse_movement || ((x!=last_x) || (y!=last_y));

  //puts("MouseRotate");
  glGetIntegerv(GL_VIEWPORT, viewp);

draw_cursor(x, y);

  if (!transversal_rotation(x, y, last_x, last_y))
    rotate((x-last_x)*sensitiveness, (y-last_y)*sensitiveness, 0.0);
  else
    {
      int s;

      if (atan2(y-(viewp[1]+viewp[3])*0.5,
		x-(viewp[0]+viewp[2])*0.5)>
	  atan2(last_y-(viewp[1]+viewp[3])*0.5,
		last_x-(viewp[0]+viewp[2])*0.5))
	s=-1;
      else
	s=1;
      
      rotate(0.0, 0.0,
	     s*sqrt((double)(SQR(x-last_x)+SQR(y-last_y)))*sensitiveness);
    }
  
  last_x0=last_x;
  last_y0=last_y;
  last_x=x;
  last_y=y;
}






void MouseZoom(int x, int y)
{
  scale(1.0+(y-last_y)*0.01, 1.0+(y-last_y)*0.01, 1.0+(y-last_y)*0.01);
  last_x=x;
  last_y=y;
draw_cursor(x, y);
}






void MouseTranslate(int x, int y)
{
  translate(xtrans+(x-last_x)*0.01, ytrans-(y-last_y)*0.01, ztrans);
  last_x=x;
  last_y=y;
draw_cursor(x, y);
}






void Mouse(int butt, int state, int x, int y)
{
  //
  // 990824: Forsoek paa aa emulere midterste tast vha. begge de to
  //         paa PC...
  //
#ifdef MIDDLE_EMU
  static int left_down=0, right_down=0;
#endif
  
/*
  printf("butt=%d state=%d x=%d y=%d "
	 "last_x=%d last_y=%d last_x=%d last_y=%d\n",
	 butt, state, x, y, last_x, last_y, last_x0, last_y0);
#ifdef MIDDLE_EMU
  printf("Left down=%d, right down=%d\n", left_down, right_down);
#endif
*/

  if (state==GLUT_DOWN)
    {
      /* A fresh start, all "last" info to be discarded... */
      last_x=x;
      last_y=y;
      last_x0=x;
      last_y0=y;

#ifdef MIDDLE_EMU
      //
      // 990824: Forsoek paa aa emulere midterste tast vha. begge de to
      //         paa PC...
      //
      if (butt==GLUT_RIGHT_BUTTON)
	{
	  if (left_down)
	    butt=GLUT_MIDDLE_BUTTON;
	  right_down=1;
	}
      else
	if (butt==GLUT_LEFT_BUTTON)
	  {
	    if (right_down)
	      butt=GLUT_MIDDLE_BUTTON;
	    left_down=1;
	  }

      /*      if (((left_down) && (butt==GLUT_RIGHT_BUTTON)) ||
	  ((right_down) && (butt==GLUT_LEFT_BUTTON)))
	{
	  butt=GLUT_MIDDLE_BUTTON;
	  left_down=right_down=1;
	}*/
#endif

      switch (butt)
	{
	case GLUT_LEFT_BUTTON:
	  xrot=0;
	  yrot=0;
	  zrot=0;
#ifndef QTMODE
	  glutMotionFunc(MouseRotate);
#else
	  // 030917: setting global variable so that new qt event handler
	  //         knows what to do.
	  extern bool enableMouseRotation;
	  enableMouseRotation=true;
#endif
	  break;
	case GLUT_MIDDLE_BUTTON:
#ifndef QTMODE
	  glutMotionFunc(MouseZoom);
#else
	  // 030918: setting global variable so that new qt event handler
	  //         knows what to do.
	  extern bool enableMouseZoom;
	  enableMouseZoom=true;
#endif
	  break;
	case GLUT_RIGHT_BUTTON:
#ifndef QTMODE
	  glutMotionFunc(MouseTranslate);
#else
	  // 030918: setting global variable so that new qt event handler
	  //         knows what to do.
	  extern bool enableMouseTranslate;
	  enableMouseTranslate=true;
#endif
	  break;
	}
    }
  
  if (state==GLUT_UP)
    {

#ifdef MIDDLE_EMU
      if (butt==GLUT_LEFT_BUTTON)
	left_down=0;
      if (butt==GLUT_RIGHT_BUTTON)
	right_down=0;
#endif
#ifdef QTMODE
      extern bool enableMouseRotation, enableMouseZoom, enableMouseTranslate;
      enableMouseRotation=false;
      enableMouseZoom=false;
      enableMouseTranslate=false;
#endif

      //
      // 020410: Introduced 'spin_speed', was 1 before, setting to 2.
      //         There seems to be some confusion here as to whether
      //         the ?rot-variables are ints or floats?!
      // 040220: Reducing to 1, seems to be very unwieldy with the new
      //         GF4 card...(?!)
      // 040505: Nå med fx5900, 1.0 ser ut til å være alt for mye for
      //         endel bview2-greier. 
      //
      //const double spin_speed=2.0;
      //const double spin_speed=1.0;
      const double spin_speed=0.1;
      
      //
      // 991123: This was a stupid thing to do! We don't need to
      //         disable this one! It's only active while the button
      //         is pressed down, anyway!
      //
      // glutMotionFunc(MouseDummy);
      if ((butt==GLUT_LEFT_BUTTON) &&
	  ((x-last_x0)*(x-last_x0)+(y-last_y0)*(y-last_y0)>4))
	{
	  if (!transversal_rotation(x, y, last_x0, last_y0))
	    {
	      xrot=(y-last_y0)*spin_speed;
	      yrot=(x-last_x0)*spin_speed;
	      /*	      zrot=0;*/
	    }
	  else
	    {
	      int s;
	      int viewp[4];

	      glGetIntegerv(GL_VIEWPORT, viewp);
	      if (atan2(y-(viewp[1]+viewp[3])*0.5,
			x-(viewp[0]+viewp[2])*0.5)>
		  atan2(last_y0-(viewp[1]+viewp[3])*0.5,
			last_x0-(viewp[0]+viewp[2])*0.5))
		s=-1;
	      else
		s=1;
	      
	      /*	      xrot=0;
	      yrot=0;*/
	      zrot=(s*(int)(sqrt((double)(SQR(y-last_y0)+SQR(x-last_x0))))*
		    spin_speed);
	    }
	}
    }	  

/*
  printf("butt=%d state=%d x=%d y=%d "
	 "last_x=%d last_y=%d last_x=%d last_y=%d\n",
	 butt, state, x, y, last_x, last_y, last_x0, last_y0);
#ifdef MIDDLE_EMU
//  printf("Left down=%d, right down=%d\n", left_down, right_down);
#endif
  puts(" ");
*/
draw_cursor(x, y);
}
