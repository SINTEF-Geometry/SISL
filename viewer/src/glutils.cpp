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

#ifdef _MSC_VER
#  include <io.h>
#else
#  include <unistd.h>
#endif

#include <fcntl.h>
#include <stdlib.h>
#include <stdio.h>

#include "glutils.h"

#ifndef M_PI
#define M_PI 3.1415926535897932384
#endif





void assert_gl_m(int line, char *file, const bool do_exit /* =true */ )
{
  int tmp=glGetError(), known_error=1;

  if (tmp==GL_NO_ERROR)
    return;
  
  while (tmp!=GL_NO_ERROR)
    {
      printf("'glGetError' returned %d ", tmp);
      switch (tmp)
	{
	case GL_NO_ERROR:
	  printf("(GL_NO_ERROR) ");
	  break;
	case GL_INVALID_ENUM:
	  printf("(GL_INVALID_ENUM) ");
	  break;
	case GL_INVALID_VALUE:
	  printf("(GL_INVALID_VALUE) ");
	  break;
	case GL_INVALID_OPERATION:
	  printf("(GL_INVALID_OPERATION) ");
	  break;
	case GL_STACK_OVERFLOW:
	  printf("(GL_STACK_OVERFLOW) ");
	  break;
	case GL_STACK_UNDERFLOW:
	  printf("(GL_STACK_UNDERFLOW) ");
	  break;
	case GL_OUT_OF_MEMORY:
	  printf("(GL_OUT_OF_MEMORY) ");
	  break;
	default:
	  known_error=0;
	  printf("(UNKNOWN ERROR(?!)) ");
	}
      if (file!=NULL)
	printf(" at line %d in file '%s'.\n", line, file);
      
      if (known_error)
	{
	  printf("The man page has this to say about this error:\n\n");
	  switch (tmp)
	    {
	    case GL_INVALID_ENUM:
	      printf("An unacceptable value is\n");
	      printf("specified for an enumerated\n");
	      printf("argument.  The offending\n");
	      printf("command is ignored, and has no\n");
	      printf("other side effect than to set\n");
	      printf("the error flag.\n");
	      break;
	    case GL_INVALID_VALUE:
	      printf("A numeric argument is out of\n");
	      printf("range.  The offending command\n");
	      printf("is ignored, and has no other\n");
	      printf("side effect than to set the\n");
	      printf("error flag.\n");
	      break;
	    case GL_INVALID_OPERATION:
	      printf("The specified operation is not\n");
	      printf("allowed in the current state.\n");
	      printf("The offending command is\n");
	      printf("ignored, and has no other side\n");
	      printf("effect than to set the error\n");
	      printf("flag.\n");
	      break;
	    case GL_STACK_OVERFLOW:
	      printf("This command would cause a\n");
	      printf("stack overflow.  The offending\n");
	      printf("command is ignored, and has no\n");
	      printf("other side effect than to set\n");
	      printf("the error flag.\n");
	      break;
	    case GL_STACK_UNDERFLOW:
	      printf("This command would cause a\n");
	      printf("stack underflow.  The\n");
	      printf("offending command is ignored,\n");
	      printf("and has no other side effect\n");
	      printf("than to set the error flag.\n");
	      break;
	    case GL_OUT_OF_MEMORY:
	      printf("There is not enough memory\n");
	      printf("left to execute the command.\n");
	      printf("The state of the GL is\n");
	      printf("undefined, except for the\n");
	      printf("state of the error flags,\n");
	      printf("after this error is recorded.\n");
	      break;
	    }
	}
      tmp=glGetError();
    }
  
  if (do_exit)
    exit(-1);
  else
    puts("\n\n\nOpenGL error!!!! \n\n\n");
}






void assert_gl_dummy_and_empty(void)
{
  return;
}






#ifndef _MSC_VER

//----------------------------------------------------------------------
//
// 040229: A small routine to list out the properties of a list of
//         GLXFBConfigs. (There *may* be properties not listed, I don't
//         remember exactly how complete this list is...)
//
//----------------------------------------------------------------------

void list_FBConfigs(GLXFBConfig *config, int nelements)
{
  int attr_list[]=
    {
      GLX_FBCONFIG_ID,
      GLX_VISUAL_ID ,
      GLX_BUFFER_SIZE,
      GLX_LEVEL 	,
      GLX_DOUBLEBUFFER 	,
      GLX_STEREO 	,
      GLX_AUX_BUFFERS 	,
      GLX_RENDER_TYPE 	,
      GLX_RED_SIZE 	,
      GLX_GREEN_SIZE 	,
      GLX_BLUE_SIZE 	,
      GLX_ALPHA_SIZE 	,
      GLX_DEPTH_SIZE 	,
      GLX_STENCIL_SIZE 	,
      GLX_ACCUM_RED_SIZE 	,
      GLX_ACCUM_GREEN_SIZE 	,
      GLX_ACCUM_BLUE_SIZE 	,
      GLX_ACCUM_ALPHA_SIZE 	,
      GLX_DRAWABLE_TYPE 	,
      GLX_X_RENDERABLE 	,
      GLX_X_VISUAL_TYPE ,
      GLX_CONFIG_CAVEAT ,
      GLX_TRANSPARENT_TYPE,
      GLX_TRANSPARENT_INDEX_VALUE 	,
      GLX_TRANSPARENT_RED_VALUE 	,
      GLX_TRANSPARENT_GREEN_VALUE 	,
      GLX_TRANSPARENT_BLUE_VALUE 	,
      GLX_TRANSPARENT_ALPHA_VALUE 	,
      GLX_MAX_PBUFFER_WIDTH 	,
      GLX_MAX_PBUFFER_HEIGHT 	,
      GLX_MAX_PBUFFER_PIXELS
    };
  char *attr_name[]=
    {
     "GLX_FBCONFIG_ID",
     "GLX_VISUAL_ID",
     "GLX_BUFFER_SIZE",
     "GLX_LEVEL",
     "GLX_DOUBLEBUFFER",
     "GLX_STEREO",
     "GLX_AUX_BUFFERS",
     "GLX_RENDER_TYPE",
     "GLX_RED_SIZE",
     "GLX_GREEN_SIZE",
     "GLX_BLUE_SIZE",
     "GLX_ALPHA_SIZE",
     "GLX_DEPTH_SIZE",
     "GLX_STENCIL_SIZE",
     "GLX_ACCUM_RED_SIZE",
     "GLX_ACCUM_GREEN_SIZE",
     "GLX_ACCUM_BLUE_SIZE",
     "GLX_ACCUM_ALPHA_SIZE",
     "GLX_DRAWABLE_TYPE",
     "GLX_X_RENDERABLE",
     "GLX_X_VISUAL_TYPE",
     "GLX_CONFIG_CAVEAT",
     "GLX_TRANSPARENT_TYPE",
     "GLX_TRANSPARENT_INDEX_VALUE",
     "GLX_TRANSPARENT_RED_VALUE",
     "GLX_TRANSPARENT_GREEN_VALUE",
     "GLX_TRANSPARENT_BLUE_VALUE",
     "GLX_TRANSPARENT_ALPHA_VALUE",
     "GLX_MAX_PBUFFER_WIDTH",
     "GLX_MAX_PBUFFER_HEIGHT",
     "GLX_MAX_PBUFFER_PIXELS"
    };
#if 0
  char *attr_str[]=
    {
      "This attribute is the XID of the GLX FBConfig.",
      "This attribute is the XID of the X Visual associated with the GLX FBConfig.",
      "This attribute defines the number of bits per color buffer. For GLX FBConfigs that correspond to a PseudoColor or StaticColor visual, this is equal to the depth value reported in the X11 visual. For GLX FBConfigs that correspond to TrueColor or DirectColor visual, this is the sum of GLX_RED_SIZE, GLX_GREEN_SIZE, GLX_BLUE_SIZE, and GLX_ALPHA_SIZE.",
      "This attribute defines the frame buffer level of the visual. Level 0 is the default frame buffer. Positive levels correspond to frame buffers that overlay the default buffer; negative levels correspond to frame buffers that underlay the default buffer.",
      "This attribute is True if color buffers exist in front/back pairs that can be swapped. Otherwise, it is False.",
      "This attribute is True if color buffers exist in left/right pairs. Otherwise, it is False.",
      "This attribute defines the number of auxiliary color buffers available. Zero indicates that no auxiliary color buffers exist.",
      "This attribute indicates what type of GLX Context a drawable created with the corresponding GLX FBConfig can be bound to. The following bit settings can exist:",
      "This attribute defines the number of red bits stored in each color buffer. If the GLX_RGBA_BIT is not set in the GLX_RENDER_TYPE attribute, the GLX_RED_SIZE attribute is undefined.",
      "This attribute defines the number of green bits stored in each color buffer. If the GLX_RGBA_BIT is not set in the GLX_RENDER_TYPE attribute, the GLX_GREEN_SIZE attribute is undefined.",
      "This attribute defines the number of blue bits stored in each color buffer. If the GLX_RGBA_BIT is not set in the GLX_RENDER_TYPE attribute, the GLX_BLUE_SIZE attribute is undefined.",
      "This attribute defines the number of alpha bits stored in each color buffer. If the GLX_RGBA_BIT is not set in the GLX_RENDER_TYPE attribute, the GLX_ALPHA_SIZE attribute is undefined.",
      "This attribute defines the number of bits in the depth buffer.",
      "This attribute defines the number of bits in the stencil buffer.",
      "This attribute defines the number of red bits stored in the accumulation buffer.",
      "This attribute defines the number of green bits stored in the accumulation buffer.",
      "This attribute defines the number of blue bits stored in the accumulation buffer.",
      "This attribute defines the number of alpha bits stored in the accumulation buffer.",
      "This attribute defines which GLX drawables are supported by the GLX FBConfig. The following bit settings can exist:",
      "This attribute indicates whether X can be used to render into a drawable created with the GLX FBConfig. This attribute is True is the GLX FBConfig supports GLX windows and/or pixmaps, otherwise it is False.",
      "This attribute defines the X visual type of the X visual associated with the GLX FBConfig. It can have one of the following values:",
      "This attribute defines any problems that the GLX FBConfig may have:",
      "This attribute defines the type of transparency (if any) supported by the FBConfig. It can have the following values:",
      "This attribute defines the index value of the transparent pixel when the transparency type is GLX_TRANSPARENT_INDEX.",
      "This attribute defines the red value of the transparent pixel when the transparency type is GLX_TRANSPARENT_RGB.",
      "This attribute defines the green value of the transparent pixel when the transparency type is GLX_TRANSPARENT_RGB",
      "This attribute defines the blue value of the transparent pixel when the transparency type is GLX_TRANSPARENT_RGB.",
      "This attribute defines the alpha value of the transparent pixel when the transparency type is GLX_TRANSPARENT_RGB.",
      "This attribute defines the maximum width value that can be passed into glXCreatePbuffer.",
      "This attribute defines the maximum height value that can be passed into glXCreatePbuffer.",
      "This attribute defines the maximum number of pixels (width times height) for a GLX Pbuffer. It can have a value that is less than the maximum width times the maximum height. Also, the value is static and assumes that no other pbuffers or X resources are contending for the framebuffer memory. Therefore, it may not be possible to allocate a pbuffer of the size given by this attribute."
    };
#endif
  int i, j;
  Display *dpy=glXGetCurrentDisplay();
  
  for (i=0; i<nelements; i++)
    {
      printf("FBConfig %d:\n", i);
      for (j=0; j<31; j++)
	{
	  int val;
	  int tmp=glXGetFBConfigAttrib(dpy, config[i], attr_list[j], &val);
	  printf("  %2d, %5d = %s: %d, tmp=%d\n",
		 j, attr_list[j], attr_name[j], val, tmp);
	}
    }

}

#endif // of if not microsoft






//======================================================================
//
// 050123: Moving this out of sfviewer1b, for re-use.
//
//         There are a set of spheres stored here, indexed by 'slot',
//         starting from 1 and onwards.
//         If 'slot' > existing number of spheres, a new one is 
//         tesselated and stored, then drawn.
//
//======================================================================

void draw_sphere(const vector3t<float> &pos, const float r,
		 const int n, const int slot,
		 const vector3t<float> &col,
		 const bool wiremode)
{
//  const int n=20; // gir 1520 triangler
//  const int n=15; // gir 840 triangler
//  const int n=29; // 3248 triangler, max som passer i grafikkortet?!
//  const int n=10; // gir 360 triangler
  //const int n=25;
//  const int n=4;
  
  static vector< vector< vector3t<float> > > all_vert1, all_norm1;
  static vector< vector<int> > all_tri;

  if (slot>(int)all_vert1.size())
    {
      vector< vector3t<float> > vert1, norm1;
      vector<int> tri;

      vert1.push_back(pos+vector3t<float>(0, 0, r));
      vert1.push_back(pos+vector3t<float>(0, 0, -r));
      norm1.push_back(vector3t<float>(0.0, 0.0, 1.0));
      norm1.push_back(vector3t<float>(0.0, 0.0, -1.0));
      
      int i, j;
      for (i=1; i<=n; i++)
	{
	  double ph=(double)i/n*M_PI-M_PI*0.5;
	  
	  for (j=0; j<2*n; j++)
	    {
	      double th=(double)j/(2*n)*M_PI*2.0;
	      
	      // sphere
	      double x=cos(th)*cos(ph);
	      double y=sin(th)*cos(ph);
	      double z=sin(ph);
	      if (i<n)
		{
		  vert1.push_back(pos+r*vector3t<float>(x, y, z));
		  norm1.push_back(vector3t<float>(x, y, z));
		}
	      
	      
	      // shall not make vertices for the top pole, but we need
	      // i=n to make the last ring of triangles!
	      
	      int k=vert1.size()-1;
	      if (j>0)
		if (i==1) {
		  tri.push_back(1); tri.push_back(k); tri.push_back(k-1);
		} else {
		  if (i==n) {
		    tri.push_back(0); tri.push_back(k-(2*n-1)+j-1);
		    tri.push_back(k-(2*n-1)+j);
		  } else {
		    tri.push_back(k-2*n-1); tri.push_back(k);
		    tri.push_back(k-1); tri.push_back(k-2*n-1);
		    tri.push_back(k-2*n); tri.push_back(k);
		  }
		}
	      else
		if (i==1) {
		  tri.push_back(1); tri.push_back(k); tri.push_back(k-1+2*n);
		} else
		  if (i==n) {
		    tri.push_back(0); tri.push_back(k);
		    tri.push_back(k-(2*n-1));
		  } else {
		    tri.push_back(k+2*n-1);	tri.push_back(k-1);
		    tri.push_back(k); tri.push_back(k);
		    tri.push_back(k-1); tri.push_back(k-2*n);
		  }
	    }
	}
      
      all_vert1.push_back(vert1);
      all_norm1.push_back(norm1);
      all_tri.push_back(tri);
    }

  if (!wiremode)
    {
      glBegin(GL_TRIANGLES);
      glColor3fv(col.raw());
      int i;
      for (i=0; i<(int)all_tri[slot-1].size(); i++)
	{
	  glNormal3fv(all_norm1[slot-1][all_tri[slot-1][i]].raw());
	  glVertex3fv(all_vert1[slot-1][all_tri[slot-1][i]].raw());
	}
      glEnd();
    }
  else
    {
      GLboolean light_temp;
      glGetBooleanv(GL_LIGHTING, &light_temp);
      glDisable(GL_LIGHTING);
      glColor3fv(col.raw());
      int i;
      for (i=0; i<(int)all_tri[slot-1].size(); i+=3)
	{
	  glBegin(GL_LINE_STRIP);
	  glVertex3fv(all_vert1[slot-1][all_tri[slot-1][i  ]].raw());
	  glVertex3fv(all_vert1[slot-1][all_tri[slot-1][i+1]].raw());
	  glVertex3fv(all_vert1[slot-1][all_tri[slot-1][i+2]].raw());
	  glVertex3fv(all_vert1[slot-1][all_tri[slot-1][i  ]].raw());
	  glEnd();
	}
      if (light_temp)
	glEnable(GL_LIGHTING);
    }
}
