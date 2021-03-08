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
#include <math.h>
#include <stdlib.h>	// atoi
#include <ctype.h>	// tolower
#include <string>
#include <iostream>
#include <fstream>
#include <stdexcept>
#include <vector>
#include <cstring>

#include <GL/gl.h>
#include <GL/glut.h>

#include "gl_aux.h"
#include "transfutils.h"
#include "mouse.h"
#include "sisl_aux.h"
//#include "read_nurbs_sf.h"
//#include "read_curves.h"
#include "GoReadWrite.h" // from 'streaming'
#include "aux2.h"
#include "jonvec.h"

#ifndef WIN32
#  include <unistd.h>	// usleep
#else
#  define MIDDLE_EMU
#endif

#ifdef _MSC_VER
#  include <time.h>
#endif

static const int DEFAULT_REF = 30000000;
static const int DEFAULT_MAX_REF = 100;

char init_key_string[1000];
double axis_thickness=0.3, axis_length=0.7, marker_size=1.0;
int draw_edges=0;
float wire_col=0.0, background_col=0.0;
int draw_axes=1;
int frames_without_movement=0;

#define MAX_SURFACES 100
#define MAX_CURVES 2000
#define MAX_LINES 100

SISLSurf *surface[MAX_SURFACES];
int surfaces=0, selected_surface=-1, surface_highlights=0, curve_highlights=0;
int selected_curve=-1;
double *normal[MAX_SURFACES];
int surf_enabled[MAX_SURFACES];
int curve_enabled[MAX_CURVES];
char *surface_name[MAX_SURFACES];
char *curve_name[MAX_CURVES];
SISLCurve *curve[MAX_CURVES];
int curves=0;
double *discr_curve[MAX_CURVES]; // This must be filled with NULLs.
vector< vector3t<float> > pcloud;

// 001101: Moving these out of set_material, so they can be used by the
//         curve-drawing functions too...
// 021203: There is a problem with duplicated code here... Should be
//         resolved some time...
//
const int predef_colours=6;
const double col_setting[3*predef_colours]
={1.0, 0.0, 0.0,
  0.0, 1.0, 0.0,
  0.0, 0.0, 1.0,
  1.0, 1.0, 0.0,
  0.0, 1.0, 1.0,
  1.0, 0.0, 1.0};
const double col_x_setting_back[3*predef_colours]
={0.6, 0.2, 0.0,
  0.0, 0.6, 0.2,
  0.2, 0.0, 0.6,
  0.6, 0.6, 0.2,
  0.2, 0.6, 0.6,
  0.6, 0.2, 0.6};


static std::string general_help_string = 
//-------|---------|---------|---------|---------|---------|---------|---------|
"USAGE:\n"
"sisl_view_demo <option list> \n"
"The options are: \n"
"'s' - next string is the filename of one or more surfaces to read \n"
"'c' - next string is the filename of one or more curves to read \n"
"'p' - next string is the filename of a pointcloud to read \n"
//"'r' - set refinement factor (Note! Must appear before 's' option.)\n"
//"'R' - set max refinement factor (Note! Must appear before 's' option.)\n"
"'r' - set max refinement factor (Note! Must appear before 's' option.)\n"
"'e' - a string with keypresses to execute follows \n"
"'hotkeys' - display a list of hotkeys that can be used \n"
"            when viewing an object\n\n";


static std::string hotkey_help_string = 
//-------|---------|---------|---------|---------|---------|---------|---------|
"\nCommand keys in viewer are: \n"
"\n"
"---General---\n"
"'q' - Quit\n"
"\n"
"---Selection---\n"
"<space> - select curve by cycling through each of them\n"
"<tab>   - select surface by cycling through each of them\n"
"\n"
"---Toggles---\n"
"'B' - toggle between black and white color for backgrounds\n"
"'A' - toggle axes on/off\n"
"'w' - toggle wireframe on/off\n"
"\n"
"---Enabling/disabling---\n"
"'e'      - enable/disable currently selected surface\n"
"<ctrl-e> - enable/disable currently selected curve\n"
"'a'      - enable all surfaces\n"
"<ctrl-a> - enable all curves\n"
"'d'      - disable all other surfaces than this\n"
"<ctrl-d> - disable all other curves than this\n"
"\n"
"---Centering of elements---\n"
"'O' - center and rescale (i.e. ditch aspect ratio) all objects around origo\n"
"      such that they fit snuggly into a (-1, -1, -1) to (1, 1, 1) cube.\n"
"'o' - center all objects around origo, no rescaling done, aspect ratios\n"
"      preserved.\n"
"\n"
"---Size of elements---\n"
"'+' - increase thickness of axes\n"
"'-' - decrease thickness of axes\n"
"'>' - increase size of points\n"
"'<' - decrease size of points\n"
"'*' / increase length of axes\n"
"'/' - decrease length of axes\n"
"\n"
"---Load/save viewpoint---\n"
"[n] is here an integer from 0 to 9.\n"
"<esc>-'w'-[n] - store viewpoint in slot [n]\n"
"<esc>-'r'-[n] - load viewpoint from slot[n]\n"
"\n\n";

static std::istream& eatwhite(std::istream& is) 
{
  char c;
  while (is.get(c)) {
    if (!isspace(c)) {
      is.putback(c);
      break;
    }
  }
  return is;
}

static void set_material(int i)
{
  const double a=0.8;
  int j;
  i=i+4;
  i=i % predef_colours;
  
  float front_mat_shininess[] = {150.0};
  float front_mat_specular[] = {0.9, 0.9, 0.9, a};

  float front_mat_ambient[] = {0.2, 0.1, 0.1, a};
  float front_mat_diffuse[] = {0.6, 0.2, 0.2, a};
  float front_mat_emission[] = {0.2, 0.05, 0.05, a};
  float back_mat_ambient[] = {0.2, 0.2, 0.1, a};
  float back_mat_diffuse[] = {0.6, 0.6, 0.2, a};
  float back_mat_emission[] = {0.2, 0.2, 0.05, a};

  for (j=0; j<3; j++)
    {
      front_mat_ambient[j] =0.1 *col_setting[     3*i+j]+0.1 ;
      front_mat_diffuse[j] =0.5 *col_setting[     3*i+j]+0.2 ;
      front_mat_emission[j]=0.15*col_setting[     3*i+j]+0.05;
      back_mat_ambient[j]  =0.1 *col_x_setting_back[3*i+j]+0.1 ;
      back_mat_diffuse[j]  =0.5 *col_x_setting_back[3*i+j]+0.2 ;
      back_mat_emission[j] =0.15*col_x_setting_back[3*i+j]+0.05;
    }
  
  glMaterialfv(GL_FRONT, GL_SHININESS, front_mat_shininess);
  glMaterialfv(GL_FRONT, GL_SPECULAR,  front_mat_specular);
  glMaterialfv(GL_FRONT, GL_DIFFUSE,   front_mat_diffuse);
  glMaterialfv(GL_FRONT, GL_AMBIENT,   front_mat_ambient);
  glMaterialfv(GL_FRONT, GL_EMISSION,  front_mat_emission);
  
  glMaterialfv(GL_BACK, GL_SHININESS, front_mat_shininess);
  glMaterialfv(GL_BACK, GL_SPECULAR,  front_mat_specular);
  glMaterialfv(GL_BACK, GL_DIFFUSE,   back_mat_diffuse);
  glMaterialfv(GL_BACK, GL_AMBIENT,   back_mat_ambient);
  glMaterialfv(GL_BACK, GL_EMISSION,  back_mat_emission);
}






static void show_general_help()
{
  puts(general_help_string.c_str());
}






static void show_hotkey_help()
{
  puts(hotkey_help_string.c_str());
}






static void draw_all_surfaces(void)
{

  int i, j, k;

  for (i=0; i<surfaces; i++)
    {
      if ((surface_highlights) && (selected_surface==i))
	glDisable(GL_LIGHTING);
      set_material(i);

      //
      // As a first brute force solution, we draw single triangles,
      // and we draw the control polygon. Normals should already be
      // computed and stored in the array 'normal'.
      //
      // A surface is drawn if it is enabled, or if it is selected.
      //
      if ((surf_enabled[i]) ||
	  ((surface_highlights) && (selected_surface==i)))
	//
	// ... or else, we fall back on the control polygon solution.
	//
	for (j=0; j<surface[i]->in1-1; j++)
	  {
	    //
	    // Now, we do a strip in the second (v)
	    // parameter direction.
	    //
	    glBegin(GL_TRIANGLE_STRIP);
	    glColor3f(1.0, 1.0, 1.0);
	    for (k=0; k<surface[i]->in2; k++)
	      {
		int p=3*(j+k*surface[i]->in1);
		glNormal3f(normal[i][p+0], normal[i][p+1], normal[i][p+2]);
		glColor3f(0.0, 1.0, 0.0);
		glTexCoord2d((double)k, 0.0);
		glVertex3f(surface[i]->ecoef[p+0],
			   surface[i]->ecoef[p+1],
			   surface[i]->ecoef[p+2]);
		p+=3;
		glNormal3f(normal[i][p+0], normal[i][p+1], normal[i][p+2]);
		glColor3f(1.0, 0.0, 0.0);
		glTexCoord2d((double)k, 1.0);
		glVertex3f(surface[i]->ecoef[p+0],
			   surface[i]->ecoef[p+1],
			   surface[i]->ecoef[p+2]);
	      }
	    glEnd();
	      
	    if (draw_edges)
	      {
		glBegin(GL_LINE_STRIP);
		glColor3f(1.0, 1.0, 1.0);
		for (k=0; k<surface[i]->in2; k++)
		  {
		    int p=3*(j+k*surface[i]->in1);
		    glNormal3f(normal[i][p+0], normal[i][p+1],
			       normal[i][p+2]);
		    glColor3f(0.0, 1.0, 0.0);
		    glTexCoord2d((double)k, 0.0);
		    glVertex3f(surface[i]->ecoef[p+0],
			       surface[i]->ecoef[p+1],
			       surface[i]->ecoef[p+2]);
		    p+=3;
		    glNormal3f(normal[i][p+0], normal[i][p+1],
			       normal[i][p+2]);
		    glColor3f(1.0, 0.0, 0.0);
		    glTexCoord2d((double)k, 1.0);
		    glVertex3f(surface[i]->ecoef[p+0],
			       surface[i]->ecoef[p+1],
			       surface[i]->ecoef[p+2]);
		  }
		  
		glEnd();
	      }

	  }
      
      if ((surface_highlights) && (selected_surface==i))
	{
	  surface_highlights--;
	  glEnable(GL_LIGHTING);
	}
      glDisable(GL_BLEND);
    }
}






static void draw_all_points(void)
{
  int i;

  glDisable(GL_LIGHTING);
  glPointSize(marker_size);
  glBegin(GL_POINTS);
  glColor3f(1.0, 1.0, 1.0);
  for (i=0; i<int(pcloud.size()); i++)
    glVertex3fv(pcloud[i].raw());
  glEnd();
  glEnable(GL_LIGHTING);
}






void draw_all_curves(void)
{
  int i, j;
  
  for (i=0; i<curves; i++)
    {
      if (curve_enabled[i])
	{
	  //
	  // As a first brute force solution, we simply evaluate the curves
	  // in 100 points and draw straight line segments.
	  // All evaluations are done once, since the curves are not evolving
	  // in any way.
	  //
#define CURVE_EVALUATIONS 500
	  if (discr_curve[i]==NULL)
	    {
	      discr_curve[i]=new double[3*CURVE_EVALUATIONS];
	      int left=0;
	      
	      for (j=0; j<CURVE_EVALUATIONS; j++)
		{
		  double t=
		    curve[i]->et[curve[i]->ik-1]+
		    (curve[i]->et[curve[i]->in]-curve[i]->et[curve[i]->ik-1])*
		    j/(CURVE_EVALUATIONS-1.0);
		  {
		    int stat;
		    
		    s1227(curve[i], 0, t, &left, discr_curve[i]+3*j, &stat);
		    if (stat!=0)
		      CRIT_ERR(printf("s1227 returned status %d.\n", stat));
		  }
		}
	    }
	  
	  //
	  // Now we draw all the line segments.
	  //
	  glLineWidth(1.0);
	  glDisable(GL_LIGHTING);
	  glBegin(GL_LINE_STRIP);
	  {
	    // int c=i % predef_colours;
	    
	    for (j=0; j<CURVE_EVALUATIONS; j++)
	      {
		//glColor3dv(col_setting+3*c);
		//glColor3f(0.0, 0.0, 0.0);
		glColor3f(1.0-background_col,
			  1.0-background_col,
			  1.0-background_col);
		glVertex3dv(discr_curve[i]+3*j);
// 		printf("%4d: %f %f %f\n", j,   discr_curve[i][3*j],
// 		       discr_curve[i][3*j+1], discr_curve[i][3*j+2]);
	      }
	  }
	  glEnd();
	  glEnable(GL_LIGHTING);
	  
	} // end of if(curve_enabled[i]) ...
    }
}






static void draw(void)
{
  /* We don't push the matrix stack, we just keep adding rotations. */
  rotate(yrot*0.01, xrot*0.01, zrot*0.01);
  
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glColorMask(GL_TRUE, GL_TRUE, GL_TRUE, GL_TRUE);
  
  draw_all_surfaces();
  draw_all_curves();
  draw_all_points();
  
  if (draw_axes)
    draw_gl_axes_old(10, axis_length,
		     axis_thickness*0.01, axis_thickness*0.04, 0.1);
  glFlush();
  glutSwapBuffers();
  glutPostRedisplay();
  
  mouse_movement=0;
}






static void Reshape(int width, int height)
{
  glViewport(0, 0, (GLint)width, (GLint)height);
  /*
    glutPostRedisplay();
  */
}






static void center_and_scale(const int rescale)
{
  int i, j;
  double minx=1e30, miny=1e30, minz=1e30;
  double maxx=-1e30, maxy=-1e30, maxz=-1e30;
  
  for (i=0; i<surfaces; i++)
    if ((surf_enabled[i]) && (surface[i]->idim==3))
      for (j=0; j<surface[i]->in1*surface[i]->in2; j++)
	{
	  minx=MIN(minx, surface[i]->ecoef[3*j+0]);
	  miny=MIN(miny, surface[i]->ecoef[3*j+1]);
	  minz=MIN(minz, surface[i]->ecoef[3*j+2]);
	  maxx=MAX(maxx, surface[i]->ecoef[3*j+0]);
	  maxy=MAX(maxy, surface[i]->ecoef[3*j+1]);
	  maxz=MAX(maxz, surface[i]->ecoef[3*j+2]);
	}
  
  for (i=0; i<curves; i++)
    if (curve_enabled[i])
      for (j=0; j<curve[i]->in; j++)
	{
	  minx=MIN(minx, curve[i]->ecoef[3*j+0]);
	  miny=MIN(miny, curve[i]->ecoef[3*j+1]);
	  minz=MIN(minz, curve[i]->ecoef[3*j+2]);
	  maxx=MAX(maxx, curve[i]->ecoef[3*j+0]);
	  maxy=MAX(maxy, curve[i]->ecoef[3*j+1]);
	  maxz=MAX(maxz, curve[i]->ecoef[3*j+2]);
	}

  for (i=0; i<int(pcloud.size()); i++)
    {
      minx=MIN(minx, pcloud[i].x());
      miny=MIN(miny, pcloud[i].y());
      minz=MIN(minz, pcloud[i].z());
      maxx=MAX(maxx, pcloud[i].x());
      maxy=MAX(maxy, pcloud[i].y());
      maxz=MAX(maxz, pcloud[i].z());
    }
  
  // printf("%f %f\n%f %f\n%f %f\n", minx, maxx, miny, maxy, minz, maxz);

  if (minx==1e30)
    {
      puts("\n\n!Oops! No surfaces/curves "
	   "enabled for 'center_and_rescale'?!\n\n");
      return;
    }
  
  double x=-0.5*(maxx+minx);
  double y=-0.5*(maxy+miny);
  double z=-0.5*(maxz+minz);
  double xs=1.0, ys=1.0, zs=1.0;
//  printf("transl: %f %f %f\n", x, y, z);
  
  if (rescale)
    {
      xs=1.0/(maxx+x);
      ys=1.0/(maxy+y);
      zs=1.0/(maxz+z);
      if (fabs(maxx-minx)<1e-10) // We're probably in a plane...
	xs=1.0;
      if (fabs(maxy-miny)<1e-10) // We're probably in a plane...
	ys=1.0;
      if (fabs(maxz-minz)<1e-10) // We're probably in a plane...
	zs=1.0;
    }

  for (i=0; i<curves; i++)
    {
      delete discr_curve[i];
      discr_curve[i]=NULL;
      for (j=0; j<curve[i]->in; j++)
	{
	  curve[i]->ecoef[3*j+0]=(curve[i]->ecoef[3*j+0]+x)*xs;
	  curve[i]->ecoef[3*j+1]=(curve[i]->ecoef[3*j+1]+y)*ys;
	  curve[i]->ecoef[3*j+2]=(curve[i]->ecoef[3*j+2]+z)*zs;
	}
    }
  
  for (i=0; i<surfaces; i++)
    {
      if (surface[i]->ikind!=1)
	printf("\nRational (?) surface probably going to be ****** up...\n");
      for (j=0; j<surface[i]->in1*surface[i]->in2; j++)
	if (surface[i]->idim==3)
	  {
	    surface[i]->ecoef[3*j+0]=(surface[i]->ecoef[3*j+0]+x)*xs;
	    surface[i]->ecoef[3*j+1]=(surface[i]->ecoef[3*j+1]+y)*ys;
	    surface[i]->ecoef[3*j+2]=(surface[i]->ecoef[3*j+2]+z)*zs;
	  }
      // if dim!=3, nothing happens... not very graceful, this, but, ...
    }

  for (i=0; i<int(pcloud.size()); i++)
    {
      pcloud[i].coo[0]=(pcloud[i].coo[0]+x)*xs;
      pcloud[i].coo[1]=(pcloud[i].coo[1]+y)*ys;
      pcloud[i].coo[2]=(pcloud[i].coo[2]+z)*zs;
    }

//  printf("xxx scale: %f %f %f\n", xscale, yscale, zscale);
  scale(1.0/xscale, 1.0/yscale, 1.0/zscale);

  if (!rescale)
    // we don't rescale the data but try to fit it into the frustum.
    {
      xs=1.0/(maxx+x);
      ys=1.0/(maxy+y);
      zs=1.0/(maxz+z);
      if (fabs(maxx-minx)<1e-10) // We're probably in a plane...
	xs=1.0;
      if (fabs(maxy-miny)<1e-10) // We're probably in a plane...
	ys=1.0;
      if (fabs(maxz-minz)<1e-10) // We're probably in a plane...
	zs=1.0;
//       printf("frustum fit: %f %f %f      %f\n",
// 	     xs, ys, zs, std::max(zs, std::max(xs, ys)));
      scale(std::max(zs, std::max(xs, ys)),
	    std::max(zs, std::max(xs, ys)),
	    std::max(zs, std::max(xs, ys)));
    }
  
}






//
// A forward decl. here, since the 'h' function only works properly
// if the 'read_curves_and_surfaces' follows after 'Key'...
//
static void read_curves_and_surfaces(int argc, char *argv[]);






static void parse_keys(const int key_in,
		       const unsigned char * const key_ptr=NULL,
		       const int keys=1)
{
  int i;
  
  for (i=0; i<keys; i++)
    {
      int rescale=0; // Used by 'o' and 'O'.
      // Default is the 'no rescale' of 'o'.
      int key;
      if (key_ptr==NULL)
	// Now 'keys' should be 1!
	key=key_in;
      else
	//
	// First, we try to detect the "pseudocharacter" '<ESC>'. Or '<TAB>'.
	// 010508: Trying to add <CTRL>
	//
	{
#define ESC_STRING "<ESC>"
#define TAB_STRING "<TAB>"
#define CTRL_STRING "<CTRL>"
	  const int esc_len=strlen(ESC_STRING);
	  const int tab_len=strlen(TAB_STRING);
	  const int ctrl_len=strlen(CTRL_STRING);

	  if ((keys-i>=esc_len+2) &&
	      (strncmp((const char *)key_ptr+i, ESC_STRING, esc_len)==0))
	    {
	      key = (27<<16) + (key_ptr[i+esc_len]<<8) + key_ptr[i+esc_len+1];
	      // 011211: Format seems to be: 27*65536+code1*256+code2.
	      i+=esc_len+1;
	    }
	  else
	    if ((keys-i>=tab_len) &&
		(strncmp((const char *)key_ptr+i, TAB_STRING, tab_len)==0))
	      {
		key=9;
		i+=tab_len-1;
	      }
	    else
	      if ((keys-i>=ctrl_len+1) &&
		  (strncmp((const char *)key_ptr+i, CTRL_STRING, ctrl_len)==0))
		{
		  key = tolower(key_ptr[i+ctrl_len])-'a'+1;
		  i+=ctrl_len;
		}
	      else
		key=key_ptr[i];
	}
      
      // printf("key=%d (%c)\n", key, key);

      switch (key)
	{
	  //
	  // 021203: What was this? Only place obj_trans (from old mouse.h or something
	  //         was used...?!)
	  //
	  //  	case 0: // @H@ Toggle mouse-translation mode.
	  //  	  obj_trans=1-obj_trans;
	  //  	  printf("obj_trans=%d\n", obj_trans);
	  //  	  break;
	  
	  //--------------------------------------------------
	  //
	  // Selection, enabling etc. of surfaces.
	  //
	  //--------------------------------------------------
	  
	case 9:
	  // (tab) Select surface by cycling through all surfaces.
	  printf("Surfaces=%d\n", surfaces);
	  if ((surfaces>0) /* && (surface_highlights==0) */ )
	    {
	      selected_surface++;
	      if (selected_surface==surfaces)
		selected_surface=0;
	      printf("Surface selected is '%s', which has dim=%d.\n",
		     surface_name[selected_surface],
		     surface[selected_surface]->idim);
	      //
	      // Now, we signal highlighting of the
	      // surface for a few frames.
	      //
	      surface_highlights=50;
	    }
	  break;

	case 'e':
	  // Enable/Disable the selected surface.
	  if (selected_surface>-1)
	    surf_enabled[selected_surface]=
	      1-surf_enabled[selected_surface];
	  if (surf_enabled[selected_surface])
	    printf("Enabled surface '%s'.\n", 
		   surface_name[selected_surface]);
	  else
	    printf("Disabled surface '%s'.\n",
		   surface_name[selected_surface]);
	  break;

	case 'd':
	  // Disable all other surfaces than this.
	  if (selected_surface>-1)
	    {
	      int i;
		      
	      for (i=0; i<surfaces; i++)
		surf_enabled[i]=(i==selected_surface);
	    }
	  printf("Enabled surface '%s'.\n",
		 surface_name[selected_surface]);
	  break;

	case 'a':
	  // Enable all surfaces.
	{
	  int i;
		  
	  for (i=0; i<surfaces; i++)
	    surf_enabled[i]=1;
	}
	puts("All surfaces enabled.");
	break;
		
	//--------------------------------------------------
	//
	// Selection, enabling etc. of curves.
	//
	//--------------------------------------------------
		
	case ' ':
	  // (space) Select curve by cycling through all of them.
	  if ((curves>0) && (curve_highlights==0))
	    {
	      selected_curve++;
	      if (selected_curve==curves)
		selected_curve=0;
	      printf("Curve selected is '%s'.\n",
		     curve_name[selected_curve]);
	      surface_highlights=5;
	    }
	  break;
		  
	case 5:
	  // (ctrl-e) Enable/Disable the selected curve.
	  if (selected_curve>-1)
	    curve_enabled[selected_curve]=1-curve_enabled[selected_curve];
	  if (curve_enabled[selected_curve])
	    printf("Enabled curve '%s'.\n", curve_name[selected_curve]);
	  else
	    printf("Disabled curve '%s'.\n", curve_name[selected_curve]);
	  break;
	case 4: // (ctrl-d) Disable all other curves than this.
	  if (selected_curve>-1)
	    {
	      int i;
	      
	      for (i=0; i<curves; i++)
		curve_enabled[i]=(i==selected_curve);
	    }
	  printf("Enabled curve '%s'.\n", curve_name[selected_curve]);
	  break;

	case 1: // (ctrl-a) Enable all curves.
	{
	  int i;
		  
	  for (i=0; i<curves; i++)
	    curve_enabled[i]=1;
	}
	puts("All curves enabled.");
	break;
	  
	//
	// 001206: Hmm... How should we handle the centering and rescaling 
	//         for inifinte straight lines? Could use the point closest
	//         to origo, perhaps???
	//
	case 'O': // Center and rescale all objects around origo.
	  rescale=1;
	  //
	  // Note how we simply continue into the code for 'o'...
	  //
	case 'o': // Center all objects around origo. Possibly rescale.
	  center_and_scale(rescale);
	  break;
	  
	case 'q': // Quit.
	  exit(0);
	  
	  //--------------------------------------------------
	  //
	  // Colour stuff.
	  //
	  //--------------------------------------------------
	  
	case 'b': // Toggle between black and white for wireframes.
	  wire_col=1.0-wire_col;
	  break;
	case 'B': // Toggle between black and white for the background.
	  background_col=1.0-background_col;
	  glClearColor(background_col, background_col, background_col, 0.0);
	  break;
	  
	  //--------------------------------------------------
	  
	case '+': // Increase thickness of axes.
	  axis_thickness*=1.1;
	  //printf("axis_thickness=%f\n", axis_thickness);
	  break;
	case '-': // Decrease thickness of axes.
	  if (axis_thickness>0.1)
	    axis_thickness*=0.9;
	  //printf("axis_thickness=%f\n", axis_thickness);
	  break;

	case '>': // Increase size of point markers.
	  marker_size*=1.1;
	  //printf("marker_size=%.5f\n", marker_size);
	  break;
	case '<': // Decrease size of point markers.
	  if (marker_size>0.00001)
	    marker_size*=0.9;
	  //printf("marker_size=%.5f\n", marker_size);
	  break;

	case '*': // Increase length of axes.
	  axis_length*=1.1;
	  //printf("axis_length=%f\n", axis_length);
	  break;
	case '/': // Decrease length of axes.
	  if (axis_length>0.1)
	    axis_length*=0.9;
	  //printf("axis_length=%f\n", axis_length);
	  break;
	  
	case 'A': // Toggle axes on/off.
	  draw_axes=1-draw_axes;
	  //printf("draw_axes=%d\n", draw_axes);
	  break;

	case 'w': // Toggle wireframe mode on/off, no hidden line removal.
	  draw_edges=1-draw_edges;
	  //printf("draw_edges=%d\n", draw_edges);
	  break;
	  
	  //--------------------------------------------------
	  //
	  // Reading and writing of viewing parameters.
	  //
	  //--------------------------------------------------
	  
	case (27<<16)+('w'<<8)+'1': // (Esc-w-<n>) Storing viewing parameters in slot <n>.
	case (27<<16)+('w'<<8)+'2':
	case (27<<16)+('w'<<8)+'3':
	case (27<<16)+('w'<<8)+'4':
	case (27<<16)+('w'<<8)+'5':
	case (27<<16)+('w'<<8)+'6':
	case (27<<16)+('w'<<8)+'7':
	case (27<<16)+('w'<<8)+'8':
	case (27<<16)+('w'<<8)+'9':
	{
	  int slot=(key&255)-'1'+1;

	  printf("Storing viewing parameters in slot %d.\n", slot);
	  char tmp[1000];
#ifdef _MSC_VER
	  // don't know which headerfile to use for VC++...
	  // Also: Using root instead of home directory...
	  sprintf(tmp, "view%d", slot);
#else
	  snprintf(tmp, 1000, "view%d", slot);
#endif
	  FILE *f=fopen(tmp, "w");
	  if (f==NULL)
	    CRIT_ERR(printf("Error opening file '%s' for writing.\n", tmp));
	  write_gl_matrices(f);
	  fprintf(f, "%.15e %.15e %.15e\n", xtrans, ytrans, ztrans);
	  fprintf(f, "%.15e %.15e %.15e\n", xscale, yscale, zscale);
	  fclose(f);
	}
	break;
	case (27<<16)+('r'<<8)+'1': // (Esc-r-<n>) Retrieving viewing parameters in slot <n>.
	case (27<<16)+('r'<<8)+'2':
	case (27<<16)+('r'<<8)+'3':
	case (27<<16)+('r'<<8)+'4':
	case (27<<16)+('r'<<8)+'5':
	case (27<<16)+('r'<<8)+'6':
	case (27<<16)+('r'<<8)+'7':
	case (27<<16)+('r'<<8)+'8':
	case (27<<16)+('r'<<8)+'9':
	{
	  int slot=(key&255)-'1'+1;

	  printf("Retrieving viewing parameters from slot %d.\n", slot);
	  char tmp[1000];
#ifdef _MSC_VER
	  // don't know which headerfile to use for VC++...
	  // Also: Using root instead of home directory...
	  sprintf(tmp, "/view%d", slot);
#else
	  snprintf(tmp, 1000, "view%d", slot);
#endif
	  FILE *f=fopen(tmp, "r");
	  if (f==NULL)
	    printf("Error opening file '%s' for reading, skipping.\n", tmp);
	  else
	    {
	      read_gl_matrices(f);
	      fscanf(f, "%lf %lf %lf\n", &xtrans, &ytrans, &ztrans);
	      fscanf(f, "%lf %lf %lf\n", &xscale, &yscale, &zscale);
	      fclose(f);
	    }
	}
	break;
	  
	default:
	  ;
	}
    }
}






static void Key(unsigned char key_in, int x, int y)
{
  static unsigned char multi_key[3];
  static int multi_key_length=0;
  int key=key_in;

  if ((key==27) && (multi_key_length==0))
    //
    // Esc as first key in sequence was pressed.
    //
    {
      multi_key[0]=key;
      multi_key_length=1;
      return;
    }

  if ((multi_key_length==1) && (multi_key[0]==27))
    //
    // Second key in esc sequence was pressed.
    //
    {
      if ((key=='r') || (key=='w'))
	{
	  multi_key_length=2;
	  multi_key[1]=key;
	  return;
	}
      //
      // This causes a two-key esc sequence to be parsed.
      //
      key+=256*multi_key[0];
      multi_key_length=0;
    }

  if (multi_key_length==2)
    //
    // Third key in a known three-key esc sequence was pressed.
    //
    {
      if ((multi_key[1]=='r') || (multi_key[1]=='w'))
	if ((key>='1') && (key<='9'))
	  printf("Storing/retrieving viewing parameters in slot '%c'.\n", key);
      //
      // This causes the three-key esc sequence to be parsed.
      //
      key=(multi_key[0]<<16) + (multi_key[1]<<8) + key;
      printf("mk=%d %d %d key=%d\n",
	     multi_key[0], multi_key[1], multi_key[2], key);
      multi_key_length=0;
    }
  
  //
  // Now key can be extended into max. 3 bytes...
  //
  parse_keys(key);

  glutPostRedisplay();
}






static void idle_func(void)
{
}






static void read_curves_and_surfaces(int argc, char *argv[])
{
  xscale=yscale=zscale=1.0;
  int ref=DEFAULT_REF; // Number of new knots between old ones.
  int maxref=DEFAULT_MAX_REF; // Maximal number of coeffs in any given direction.
  // (n, n) makes sure new knots are inserted close to max limit = n...
  int i;

  //
  // This must be reset every time we change control vertices, since
  // the discretization is done in the drawing routine.
  //
  for (i=0; i<curves; i++) {
    if (discr_curve[i]!=NULL) {
      delete discr_curve[i];
      discr_curve[i]=NULL;
    }
  }

  // @HH@
  // @HH@ Use the following optional command line options:
  // @HH@
  
  surfaces=0;
  curves=0;



  for (i=1; i<argc; i++) {
    switch (argv[i][0]) {
    case 's': {
      // Next string is filename for surface. (One surface.)
	    
      puts("Reading surfaces");
      // surface[surfaces-1]=read_nurbs_sf(argv[i+1]); // old format    

      std::vector<SISLSurf*> tmp;
      std::ifstream is(argv[i+1]);
      if (!is) {
	CRIT_ERR(printf("Could not open file: '%s'.\n", argv[i+1]));
      }
      try {
	eatwhite(is);
	while (!is.eof()) {
	  //surface[surfaces-1] = readGoSurface(is);
	  tmp.push_back(readGoSurface(is));
	  eatwhite(is);
	}
      } catch (std::exception& e) {
	CRIT_ERR(printf("Error occured while reading surface: %s\n", e.what()));
      }
      is.close();
      int num_surfaces = tmp.size();
      if (surfaces + num_surfaces > MAX_SURFACES) {
	CRIT_ERR(puts("Increase MAX_SURFACES."));
      }
	    
      for (int k = 0; k < num_surfaces; ++k) {
	//
	// 010116: This should be a quick fix for periodic
	//         surfaces...
	//
	if (tmp[k] == NULL) {
	  CRIT_ERR(printf("Couldn't read SISLSurf '%s'.\n", argv[i+1]));	    
	}

	if ((tmp[k]->cuopen_1 == SISL_CRV_PERIODIC ||
	     tmp[k]->cuopen_2 == SISL_CRV_PERIODIC)) {
	  int kstat;
	  SISLSurf *tmp_surf;
	  make_sf_kreg(tmp[k], &tmp_surf, &kstat);
	  if (kstat < 0) {
	    CRIT_ERR(printf("make_sf_kreg failed!\n"));
	  }
	  freeSurf(tmp[k]);
	  tmp[k] = tmp_surf;
	}
	if (tmp[k]->idim != 3) {
	  CRIT_ERR(printf("Dimension of surface is %d and not 3!\n",
			  tmp[k]->idim));
	}
		
	if (surface[surfaces]) {
	  // deleting old surface
	  freeSurf(surface[surfaces]);
	}
	surface[surfaces] = tmp[k];
	surface_name[surfaces] = argv[i+1];
	surf_enabled[surfaces] = 1;

	// Generating an approximating polygon.
	lower_degree_and_subdivide(surface + surfaces, ref, maxref);

	// evaluating normals (normalized)
	delete normal[surfaces];
	compute_surface_normals(surface[surfaces], normal + surfaces);
		
	++surfaces;
      }
    }
      i++;
      break;
      
    case 'c': {
      // Next string is filename for file containing 1 curve.
	    
      printf("Reading a single curve into slot %d.\n", curves);

      std::vector<SISLCurve*> tmp;
	    
      //n=get_curve_set(argv[i+1], &tmp, &stat);
      //get_single_curve(argv[i+1], &tmp, &stat); // old format
      std::ifstream is(argv[i+1]);
      if (!is) {
	CRIT_ERR(printf("Could not open file: '%s'.\n", argv[i+1]));
      }
      try {
	eatwhite(is);
	while (!is.eof()) {
	  tmp.push_back(readGoCurve(is));
	  eatwhite(is);
	}
      } catch (std::exception& e) {
	CRIT_ERR(printf("Error occured while reading curve: %s\n", e.what()));
      }
      is.close();
      int num_curves = tmp.size();

      if (curves + num_curves > MAX_CURVES) {
	CRIT_ERR(puts("Increase MAX_CURVES."));
      }
	    
      for(int k = 0; k < num_curves; ++k) {
	if (curve[curves + k] != NULL) {
	  freeCurve(curve[curves + k]);
	}
	curve[curves + k] = tmp[k];
		
	//
	// 001206: If the dimension is 2, we set up a
	//         new curve, filling in zeros.
	//
	if (curve[curves + k]->idim==2) {
		    
	  double *new_coeffs=
	    new double[3*curve[curves + k]->in];
	  if (new_coeffs==NULL)
	    CRIT_ERR(puts("Couldn't allocate memory."));
	  int j;
		    
	  for (j=0; j<curve[curves + k]->in; j++) {
	    new_coeffs[3*j+0]= curve[curves + k]->ecoef[2*j+0];
	    new_coeffs[3*j+1]= curve[curves + k]->ecoef[2*j+1];
	    new_coeffs[3*j+2]=0.0;
	  }
	  SISLCurve *tmp2=curve[curves + k];
	  curve[curves + k]= newCurve(tmp2->in, tmp2->ik, tmp2->et,
				      new_coeffs, tmp2->ikind, 3, 1);
	  freeCurve(tmp2);
	}
		
	if (curve[curves + k]->idim!=3) {
	  CRIT_ERR(printf("Dimension of curve is %d and not 3?!\n",
			  curve[curves + k]->idim));
	}

	curve_name[curves + k]=argv[i+1];
	curve_enabled[curves + k]=1;
      }
      curves+=num_curves;
    }
      i++;
      break;
	    
    case 'p': {
      // Next string is filename for file containing a point cloud.
	    
      printf("Reading a point cloud.\n");

      std::ifstream is(argv[i+1]);
      if (!is) {
	CRIT_ERR(printf("Could not open file: '%s'.\n", argv[i+1]));
      }
      try {
	  vector<double> coords;
	  readGoPoints(coords, is);
	  int num_points = int(coords.size()) / 3;
	  printf("Number of vertices: %d\n", num_points);
	  for (int i = 0; i < num_points; ++i) {
	      pcloud.push_back(vector3t<float>(coords[3 * i], 
					       coords[3 * i + 1], 
					       coords[3 * i + 2]));
	  }
// 	eatwhite(is);
// 	int tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	is >> tmp;
// 	eatwhite(is);
// 	is >> tmp;
// 	printf("Number of vertices: %d\n", tmp);
// 	while (!is.eof()) {
// 	  double x, y, z;
// 	  is >> x;
// 	  is >> y;
// 	  is >> z;
// 	  //printf("point: %f %f %f\n", x, y, z);
// 	  pcloud.push_back(vector3t<float>(x, y, z));
// 	  eatwhite(is);
//      }
      } catch (std::exception& e) {
	  CRIT_ERR(printf("Error occured while reading curve: %s\n", e.what()));
      }
      is.close();
      printf("pcloud size is now %d\n", pcloud.size());
    }
	i++;
      break;
	    
// 	case 'r':
// 	    // Set refinement factor. Default value is ???.
	    
// 	    puts("Reading surface refinement factor");
// 	    ref=atoi(argv[i+1]);
// 	    i++;
// 	    break;
	    
//	case 'R': 
    case 'r':
      // Set max refinement factor. Default value is ???.
	    
      puts("Reading upper bound for surface refinement.");
      maxref=atoi(argv[i+1]);
      i++;
      break;
	    
    case 'e': 
      // String with keypresses to execute follows. Not
      // everything will work!
	    
      printf("Executing '%s'.\n", argv[i+1]);
      strncpy(init_key_string, argv[i+1], 1000);
      i++;
      break;
	    
    default:
      puts("Huh?! Unknown option.");
      exit(0);
	    
    }
  }
}





int main(int argc, char *argv[])
{
  if (argc == 1) {
    // program called with no arguments.  Show info
    show_general_help();
    return 0;
  } else if (argc == 2 && 
	     (std::string(argv[1]) == "hotkeys" || std::string(argv[1]) == "HOTKEYS")) {
    show_hotkey_help();
    return 0;
  }
  xscale=yscale=zscale=1.0;
  int i;

  
  for (i=0; i<MAX_SURFACES; i++)
    {
      surface[i]=NULL;
      normal[i]=NULL;
    }
  for (i=0; i<MAX_CURVES; i++)
    curve[i]=NULL;

  for (i=0; i<MAX_CURVES; i++)
    discr_curve[i]=NULL;

  read_curves_and_surfaces(argc, argv);
  
  {
    int tx=32, ty=32;	// gl_init needs variables, but doesn't do
    // anything when texfile==NULL.
#ifdef _MSC_VER
    gl_init(argc, argv, 500, 500, 180, 100, true,
#else
	    gl_init(argc, argv, 1100, 1100, 460, 20, true,
#endif
		    1, GL_TRUE, 1, GL_TRUE,
		    xtrans, ytrans, ztrans,
		    // xscale, yscale, zscale, tx, ty, NULL);
		    xscale, yscale, zscale, tx, ty,
		    // GL_MODULATE,
		    // GL_BLEND,
		    GL_DECAL,
		    NULL			// NULL produces a frame.
		    // "greyclouds3.pgm"	// A texture.
	      );
	    }
  
    //
    // This will have the same effect as pressing 'tab'+'d'.
    //
    if (surfaces>0)
      {
	int i;
      
	selected_surface=0;
	for (i=0; i<surfaces; i++)
	  surf_enabled[i]= 1 ; //(i==selected_surface);
      }

    //
    // This will have the same effect as pressing 'space'+'ctrl-d'.
    //
    if (curves>0)
      {
	int i;
      
	selected_curve=0;
	for (i=0; i<curves; i++)
	  curve_enabled[i]=1; // (i==selected_curve);
      }

    glPolygonOffset(1.0, 1.0);
    glEnable(GL_POLYGON_OFFSET_POINT);
    glEnable(GL_POLYGON_OFFSET_LINE);
    glEnable(GL_POLYGON_OFFSET_FILL);
  
    glutReshapeFunc(Reshape);
    glutKeyboardFunc(Key);
    glutMotionFunc(MouseRotate);
    glutMouseFunc(Mouse);
    glutDisplayFunc(draw);
    glutIdleFunc(idle_func);

    scale(0.5, 0.5, 0.5);
    glEnable(GL_NORMALIZE);


    parse_keys(0,
	       (const unsigned char *)init_key_string,
	       strlen(init_key_string));

    // making sure loaded object is centered
    //parse_keys('o'); // doesn't work as expected


    glutMainLoop();
  
    return 0;
  }
