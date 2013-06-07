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

#ifndef GL_AUX_INCLUDED
#define GL_AUX_INCLUDED

#include <stdio.h>
#include <GL/glut.h>
#include <algorithm>

#include "transfutils.h"
#include "jonvec.h"

extern const int predefined_colours;
extern const double predef_col_table[];
extern const double col_setting_back_old[];
extern const double col_setting_back[];


//
// 020723: I don't remember when the stuff above is used... May be that the
//         class below could/should replace them...
//
// This class is used in this way: The glut application should make a
// (global) instance of this class, and call it's 'Key' member from
// the 'Key' callback function. The "multi purpose standard
// Key-callback" defined in xxxxx does just this.
//
// 020822: This crap should be fixed, updated, cleaned up.
//

class material_appearance
{
private:
#ifndef MICROSOFT
  static const int schemes=7;
#else
  //
  // @@@ 021027: Remember to ask afr about this...
  //
#  define schemes 7
#endif

  int cubemap_mode, col_scheme;
  float front_mat_shininess[schemes];
  float front_mat_specular[4*schemes];
  float front_mat_diffuse[4*schemes];
  float front_mat_ambient[4*schemes];
  float front_mat_emission[4*schemes];
  float back_mat_shininess[schemes];
  float back_mat_specular[4*schemes];
  float back_mat_diffuse[4*schemes];
  float back_mat_ambient[4*schemes];
  float back_mat_emission[4*schemes];
  int spec_mode;
  int col_face;		// State variable for chosen face to modify.
  int col_param;	// State variable for chosen param to modify.
  float *property;	// Property data to modify.
public:
  double ztrans_delta;
private:
  int rescale, bg_invert;
  
public:
  material_appearance(void):
    cubemap_mode(0), col_scheme(0), spec_mode(0), col_face(0), col_param(0),
    property(front_mat_specular), ztrans_delta(0.0), rescale(0), bg_invert(0)
    {
      static const float front_mat_shininess0[schemes] ={200.0,
							 200.0,
							 200.0,
							 200.0,
							 200.0,
							 200.0,
							 200.0};
      static const float front_mat_specular0[4*schemes]={0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0,
							 0.8, 0.8, 0.8, 1.0};
      static const float front_mat_diffuse0[4*schemes] ={0.7, 0.7, 0.7, 1.0,
							 0.7, 0.4, 0.4, 1.0,
							 0.4, 0.7, 0.7, 1.0,
							 0.0, 0.8, 0.2, 1.0,
							 0.7, 0.7, 0.4, 1.0,
							 0.4, 0.4, 0.7, 1.0,
							 0.5, 0.4, 0.4, 1.0};
      static const float front_mat_ambient0[4*schemes] ={0.3, 0.3, 0.3, 1.0,
							 0.8, 0.8, 0.3, 1.0,
							 0.3, 0.3, 0.8, 1.0,
							 0.5, 0.0, 0.0, 1.0,
							 0.8, 0.3, 0.8, 1.0,
							 0.3, 0.8, 0.3, 1.0,
							 0.0, 0.0, 0.0, 1.0};
      static const float front_mat_emission0[4*schemes]={0.3, 0.0, 0.0, 1.0,
							 0.3, 0.2, 0.0, 1.0,
							 0.0, 0.2, 0.0, 1.0,
							 0.3, 0.2, 0.3, 1.0,
							 0.1, 0.2, 0.1, 1.0,
							 0.0, 0.2, 0.4, 1.0,
							 0.0, 0.0, 0.0, 1.0};
      static const float back_mat_shininess0[schemes]  ={200.0, 
							 200.0,
							 228.0,
							 228.0,
							 228.0,
							 228.0,
							 228.0};
      static const float back_mat_specular0[4*schemes] ={0.8, 0.0, 0.0, 1.0,
							 0.0, 0.0, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0};
      static const float back_mat_diffuse0[4*schemes]  ={0.7, 0.7, 0.7, 1.0,
							 0.6, 0.9, 0.6, 1.0,
							 0.0, 0.7, 0.9, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.9, 0.7, 0.9, 1.0,
							 0.4, 0.4, 0.4, 1.0};
      static const float back_mat_ambient0[4*schemes]  ={0.3, 0.3, 0.3, 1.0,
							 0.9, 0.3, 0.3, 1.0,
							 0.3, 0.9, 0.3, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.9, 0.9, 0.3, 1.0,
							 0.0, 0.0, 0.0, 1.0};
      static const float back_mat_emission0[4*schemes] ={0.0, 0.3, 0.0, 1.0,
							 0.9, 0.3, 0.0, 1.0,
							 0.0, 0.9, 0.9, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.8, 0.9, 0.0, 1.0,
							 0.9, 0.0, 0.9, 1.0,
							 0.0, 0.0, 0.0, 1.0};
      
      memcpy(front_mat_shininess, front_mat_shininess0, schemes*sizeof(float));
      memcpy(front_mat_specular, front_mat_specular0, 4*schemes*sizeof(float));
      memcpy(front_mat_diffuse, front_mat_diffuse0, 4*schemes*sizeof(float));
      memcpy(front_mat_ambient, front_mat_ambient0, 4*schemes*sizeof(float));
      memcpy(front_mat_emission, front_mat_emission0, 4*schemes*sizeof(float));
      memcpy(back_mat_shininess, back_mat_shininess0, 4*schemes*sizeof(float));
      memcpy(back_mat_specular, back_mat_specular0, 4*schemes*sizeof(float));
      memcpy(back_mat_diffuse, back_mat_diffuse0, 4*schemes*sizeof(float));
      memcpy(back_mat_ambient, back_mat_ambient0, 4*schemes*sizeof(float));
      memcpy(back_mat_emission, back_mat_emission0, 4*schemes*sizeof(float));
    }

  int col_scheme_selected(void) const { return col_scheme; }
  void col_scheme_select(const int i)
    {
      col_scheme=i%schemes;
    }

  void set_material(void)
    {
      glMaterialfv(GL_FRONT, GL_SHININESS, front_mat_shininess+col_scheme);
      glMaterialfv(GL_FRONT, GL_SPECULAR, front_mat_specular+col_scheme*4);
      glMaterialfv(GL_FRONT, GL_DIFFUSE, front_mat_diffuse+col_scheme*4);
      glMaterialfv(GL_FRONT, GL_AMBIENT, front_mat_ambient+col_scheme*4);
      glMaterialfv(GL_FRONT, GL_EMISSION, front_mat_emission+col_scheme*4);
      
      glMaterialfv(GL_BACK, GL_SHININESS, back_mat_shininess+col_scheme);
      glMaterialfv(GL_BACK, GL_SPECULAR, back_mat_specular+col_scheme*4);
      glMaterialfv(GL_BACK, GL_DIFFUSE, back_mat_diffuse+col_scheme*4);
      glMaterialfv(GL_BACK, GL_AMBIENT, back_mat_ambient+col_scheme*4);
      glMaterialfv(GL_BACK, GL_EMISSION, back_mat_emission+col_scheme*4);
    }

  void redefine_material_f(const vector3t<double> &diff,
			   const vector3t<double> &amb,
			   const vector3t<double> &em,
			   const vector3t<double> &spec,
			   const double shin)
    {
      front_mat_diffuse[col_scheme*4+0]=diff.x();
      front_mat_diffuse[col_scheme*4+1]=diff.y();
      front_mat_diffuse[col_scheme*4+2]=diff.z();
      front_mat_diffuse[col_scheme*4+3]=1.0;

      front_mat_ambient[col_scheme*4+0]=amb.x();
      front_mat_ambient[col_scheme*4+1]=amb.y();
      front_mat_ambient[col_scheme*4+2]=amb.z();
      front_mat_ambient[col_scheme*4+3]=1.0;

      front_mat_emission[col_scheme*4+0]=em.x();
      front_mat_emission[col_scheme*4+1]=em.y();
      front_mat_emission[col_scheme*4+2]=em.z();
      front_mat_emission[col_scheme*4+3]=1.0;

      front_mat_specular[col_scheme*4+0]=spec.x();
      front_mat_specular[col_scheme*4+1]=spec.y();
      front_mat_specular[col_scheme*4+2]=spec.z();
      front_mat_specular[col_scheme*4+3]=1.0;

      front_mat_shininess[col_scheme]=shin;
    }

  void Key(unsigned char key)
    {
      bool adjust_material=false;
      int comp=0; // Colour component to modify.
      double delta=0.0; // Amount of modification.
      
      switch (key)
	{
	  // case 'm': // @H@ ------ Material and lighting properties prefix:
	  
	case 's': // @H@ Toggle on/off specular property of material.
	  spec_mode=1-spec_mode;
	  printf("spec_mode=%d\n", spec_mode);
	  {
	    int i;
	    
	    for (i=0; i<3; i++)
	      if (spec_mode)
		front_mat_specular[4*schemes+i]=
		  back_mat_specular[4*schemes+i]=1.0;
	      else
		front_mat_specular[4*schemes+i]=
		  back_mat_specular[4*schemes+i]=0.0;
	  }
	  glMaterialfv(GL_FRONT, GL_SPECULAR,
		       front_mat_specular+4*col_scheme);
	  glMaterialfv(GL_BACK, GL_SPECULAR,
		       back_mat_specular+4*col_scheme);
	  //material_prefix_pressed=false;
	  break;
	  
	case 'c': // @H@ Cycle through color schemes.
	  col_scheme=(col_scheme+1)%schemes;
	  printf("col_scheme=%d\n", col_scheme);
	  adjust_material=true;
	  //material_prefix_pressed=false;
	  break;

	case 'b': // @H@ Toggle which face to adjust colour for.
	  col_face=1-col_face;
	  printf("col_face=%s\n", col_face ? "back" : "front");
	  break;
	  
	case 'p': // @H@ Cycle through colour parameter to tune.
	  col_param=(col_param+1)%5;
	  printf("Colour parameter to be tuned: '");
	  switch (col_param)
	    {
	    case 0:
	      printf("shininess");
	      property=col_face ? back_mat_shininess : front_mat_shininess;
	      break;
	    case 1:
	      printf("specular");
	      property=col_face ? back_mat_specular : front_mat_specular;
	      break;
	    case 2:
	      printf("diffuse");
	      property=col_face ? back_mat_diffuse : front_mat_diffuse;
	      break;
	    case 3:
	      printf("ambient");
	      property=col_face ? back_mat_ambient : front_mat_ambient;
	      break;
	    case 4:
	      printf("emission");
	      property=col_face ? back_mat_emission : front_mat_emission;
	      break;
	    }
	  printf("'.\n");
	  break;

	case 'r': // @H@ Increase red-comp. of current colour.
	  comp=0; adjust_material=true; delta=1.0/16.0; break;
	case 'f': // @H@ Decrease red-comp. of current colour.
	  comp=0; adjust_material=true; delta=-1.0/16.0; break;
	case 't': // @H@ Increase green-comp. of current colour.
	  comp=1; adjust_material=true; delta=1.0/16.0; break;
	case 'g': // @H@ Decrease green-comp. of current colour.
	  comp=1; adjust_material=true; delta=-1.0/16.0; break;
	case 'y': // @H@ Increase blue-comp. of current colour.
	  comp=2; adjust_material=true; delta=1.0/16.0; break;
	case 'h': // @H@ Decrease blue-comp. of current colour.
	  comp=2; adjust_material=true; delta=-1.0/16.0; break;

	case 'H': // @H@ Help for the 'material property' mode.


	  puts("\nUse the following keys:\n");
// 	  system("cat "PWD"/"__FILE__" | "
// 		 "gawk '/[c]ase.*@H@/ { match($0, "
// 		 "'\"/[c]ase.*'.*'/\"'); printf(\"%c) \", substr($0, "
// 		 "RSTART+RLENGTH-2, 1)); match($0, /@H@/); print(substr($0, "
// 		 "RSTART+RLENGTH+1)) } /@[H]H@/ { match($0, /@[H]H@/); "
// 		 "print(substr($0, RSTART+RLENGTH+1)) }'");


	  puts("\nSpecific help for material property selection:\n"
	       "  Note that the following keys do not exit 'material' mode\n"
	       "  after their application. (The key 'm' does that.) Use:\n"
	       "  b - for toggling of face, 'back'/'front',\n"
	       "  p - for cycling through allowed properties,\n"
	       "      'shininess', 'specular', 'diffuse', 'ambient' and\n"
	       "      'emission',\n"
	       "  r - for increasing the red-component,\n"
	       "  f - for decreasing the red-component,\n"
	       "  t - for increasing the green-component,\n"
	       "  g - for decreasing the green-component,\n"
	       "  y - for increasing the blue-component and\n"
	       "  h - for decreasing the blue-component.\n"
	       "  Also, these are available:\n"
	       "  i - invert background,"
	       "  d - dump numbers.\n");
	  break;

	case 'i': // @H@ Invert background colour.
	  bg_invert=1-bg_invert;
	  printf("bg_invert=%d\n", bg_invert);
	  if (bg_invert)
	    glClearColor(1.0, 1.0, 1.0, 0.0);
	  else
	    glClearColor(0.0, 0.0, 0.0, 0.0);
	  break;
	  
	case 'd':
	{
	  int i;
	  
	  puts("\nThe numbers for front properties, listed in the order of\n"
	       "shininess, specular, diffuse, ambient and emmisive light,\n"
	       "is:\n");
	  printf("%.5f\n", *(front_mat_shininess+col_scheme));
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (front_mat_specular+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (front_mat_diffuse+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (front_mat_ambient+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (front_mat_emission+col_scheme*4)[i]);
	  puts("\n\nThe numbers for back properties, listed in the order of\n"
	       "shininess, specular, diffuse, ambient and emmisive light,\n"
	       "is:\n");
	  printf("%.5f\n", *(back_mat_shininess+col_scheme));
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (back_mat_specular+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (back_mat_diffuse+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (back_mat_ambient+col_scheme*4)[i]);
	  printf("\n");
	  for (i=0; i<4; i++)
	    printf("%.5f   ", (back_mat_emission+col_scheme*4)[i]);
	  printf("\n\n");
	}
	break;

	case 'z': // @H@ Reduce the amount of 'perspective'.
	  if (ztrans_delta==0.0)
	    ztrans_delta=0.01;
	  ztrans_delta=std::min(ztrans_delta*1.1, 1000.0);
	  ztrans=-6-ztrans_delta;
	  printf("ztrans_delta=%f\n", ztrans_delta);
	  printf("ztrans=%f\n", ztrans);
	  glMatrixMode(GL_PROJECTION);
	  glLoadIdentity();
	  glFrustum(-1.0, 1.0, -1.0, 1.0, 5+ztrans_delta, 20+ztrans_delta);
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();
	  printf("trans=%f %f %f\n", xtrans, ytrans, ztrans);
	  glTranslated(xtrans, ytrans, ztrans);
	  printf("scale=%f %f %f\n", xscale, yscale, zscale);
	  glScaled(xscale, yscale, zscale);
	  break;
	case 'x': // @H@ Increase the amount of 'perspective'.
	  ztrans_delta*=0.9;
	  ztrans=-6-ztrans_delta;
	  printf("ztrans_delta=%f\n", ztrans_delta);
	  printf("ztrans=%f\n", ztrans);
	  glMatrixMode(GL_PROJECTION);
	  glLoadIdentity();
	  glFrustum(-1.0, 1.0, -1.0, 1.0, 5+ztrans_delta, 20+ztrans_delta);
	  glMatrixMode(GL_MODELVIEW);
	  glLoadIdentity();
	  glTranslated(xtrans, ytrans, ztrans);
	  glScaled(xscale, yscale, zscale);
	  break;
	
	case 'm': // @H@ ------ Material prefix section ends here.
	  //material_prefix_pressed=false;
	  break;
	  
	case 'q': // @H@ Exit the program.
	  exit(0);
	  break;

	}
      
      if (adjust_material)
	{
	  float * const &tmp=(property+
			      col_scheme*(col_param==0 ? 1 : 4)+
			      (col_param==0 ? 0 : comp));
	  printf("Old value %.4f replaced by new value ", *tmp);
	  if (col_param==0)
	    *tmp=std::min(std::max((*tmp)*(1+delta), 0.0), 1000.0);
	  else
	    *tmp=std::min(std::max(*tmp+delta, 0.0), 1.0);
	  printf("%.4f\n", *tmp);

	  set_material();
	}

    } // end of 'Key'

};






void draw_cylinder(double x0, double y0, double z0,
		   double x1, double y1, double z1,
		   double radius, double radius2, int n);

void draw_gl_axes_old(int n=10, double r=0.7, double radius=0.01,
		      double rim=0.04, double l=0.1);

void draw_grid(const int n1, const int n2, const int n3);

void draw_grid_planes(const int n1, const int n2, const int n3);

void gl_init(int argc, char *argv[],
	     const int xs, const int ys,
	     const int x0, const int y0,
	     const int doubleBuffer,
	     const int two_sided,
	     const int lighting,
	     const int normalize,
	     const int smooth,
	     const double xtrans, const double ytrans, const double ztrans,
	     const double xscale, const double yscale, const double zscale,
	     int &tx,
	     int &ty,
	     const int texture_mode,
	     const char * const texfile=NULL);

void print_gl_matrix(const int m);

void reshape_window(int width, int height);

//
// Writing and reading the GL transformation matrices to and from
// files.
// 020505: Copied from the sisl_view version of gl_aux... Here is some
//         more merging to do...
//
void write_gl_matrices(FILE *f);
void read_gl_matrices(FILE *f);






#endif
