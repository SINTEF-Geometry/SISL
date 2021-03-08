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

#ifndef JONVEC_H_INCLUDED
#define JONVEC_H_INCLUDED

#include <stdio.h>
#include <math.h>
#include <vector>

#include "aux2.h"


/*!

  \page filters_mainpage The filter routines main page
  \anchor filters_mainpage
  
  \section filters_intro Introduction

  This is a collection of routines for handling various (3x3x3) 3d
  filters.<p>

  020114: This REALLY REALLY must be cleaned up...
  021231: Not yet, but soon, it will only be a question of removing commented
          out stuff... Hopefully...
  050120: Changed return types for boolean functions from int to bool.
          Hope this doesn't break anything.

*/






template<typename T>
class vector3t
{
public:
  T coo[3];

  vector3t()
    {
      coo[0]=0.0;
      coo[1]=0.0;
      coo[2]=0.0;
    };
  vector3t(const T a, const T b, const T c)
    {
      coo[0]=a;
      coo[1]=b;
      coo[2]=c;
    }
  vector3t(const float *a)
    {
      coo[0]=a[0];
      coo[1]=a[1];
      coo[2]=a[2];
    }
  vector3t(const double *a)
    {
      coo[0]=a[0];
      coo[1]=a[1];
      coo[2]=a[2];
    }
  // 030819: Disse to ødelegger et eller annet. Fordi de erstatter en default (og bitwise) copy constructor, uten at de samtidig er gode nok til å gjøre det??? Skjønner ikke helt hva som skjer... Må sjekkes opp... @@@
//    vector3t(vector3t<float> &v0) // 030819
//      {
//        coo[0]=v0.x();
//        coo[1]=v0.y();
//        coo[2]=v0.z();
//      }
//    vector3t(vector3t<double> &v0) // 030819
//      {
//        coo[0]=v0.x();
//        coo[1]=v0.y();
//        coo[2]=v0.z();
//      }
  ~vector3t(void) {}

//  vector rotate2d(T,int);
//  vector transform(coo_system &);

  //
  // 010122: For fast access used in OpenGL:
  //
  inline const T *raw(void) const
    {
      return coo;
    }

  //
  // 010125
  // 021024: Added the reference-versions.
  //
  inline T x(void) const { return coo[0]; }
  inline T y(void) const { return coo[1]; }
  inline T z(void) const { return coo[2]; }
  inline T &x(void) { return coo[0]; }
  inline T &y(void) { return coo[1]; }
  inline T &z(void) { return coo[2]; }

  //
  // 010207
  //
  inline void setx(const T a) { coo[0]=a; }
  inline void sety(const T a) { coo[1]=a; }
  inline void setz(const T a) { coo[2]=a; }

  //
  // Check out this on gcc... Why does only SGI CC complain?
  //
  inline bool operator==(const vector3t &v) const
    {
      return ((coo[0]==v.coo[0]) && (coo[1]==v.coo[1]) && (coo[2]==v.coo[2]));
    }

  //
  // 050120: Adding this... Would also be nice with one "and"'ing the
  //         results on each component... No, not necessary...
  //         Remember that deMorgan gives that a<&b <=> !(a>|b).
  //
  inline bool operator<(const vector3t &v) const
    {
      return ((coo[0]<v.coo[0]) || (coo[1]<v.coo[1]) || (coo[2]<v.coo[2]));
    }
  inline bool operator>(const vector3t &v) const
    {
      return ((coo[0]>v.coo[0]) || (coo[1]>v.coo[1]) || (coo[2]>v.coo[2]));
    }

  //
  // Note: Prefer a!=b to !(a==b) since the former is probably faster!
  //
  inline bool operator!=(const vector3t &v) const
    {
      return ((coo[0]!=v.coo[0]) || (coo[1]!=v.coo[1]) || (coo[2]!=v.coo[2]));
    }

  //
  // Ok, the strange warning produced when actually using this operator,
  // stemmed from a missing 'const' at the end of the declaration of
  // this operator. (Then the function could have altered *this, and
  // that would have discarded any const'ness...)
  //
  inline const vector3t operator+(const vector3t &v) const
    {
      return vector3t(coo[0]+v.coo[0], coo[1]+v.coo[1], coo[2]+v.coo[2]);
    }

  inline vector3t &operator+=(const vector3t &v)
    {
      coo[0]+=v.coo[0]; coo[1]+=v.coo[1]; coo[2]+=v.coo[2];
      return *this;
    }

  // 030819
//    inline vector3t<float> &operator+=(const vector3t<double> &v)
//      {
//        coo[0]+=v.coo[0]; coo[1]+=v.coo[1]; coo[2]+=v.coo[2];
//        return *this;
//      }

  inline vector3t &operator-=(const vector3t &v)
    {
      coo[0]-=v.coo[0]; coo[1]-=v.coo[1]; coo[2]-=v.coo[2];
      return *this;
    }

  // 030513
  inline vector3t &operator-=(const T &d)
    {
      coo[0]-=d; coo[1]-=d; coo[2]-=d;
      return *this;
    }

  // 030626
  inline vector3t &operator+=(const T &d)
    {
      coo[0]+=d; coo[1]+=d; coo[2]+=d;
      return *this;
    }

  // 030709
  inline vector3t &operator*=(const double x)
    {
      coo[0]*=x; coo[1]*=x; coo[2]*=x;
      return *this;
    }

  // 030709
  inline vector3t &operator/=(const double x)
    {
      coo[0]/=x; coo[1]/=x; coo[2]/=x;
      return *this;
    }

  inline const vector3t operator-(const vector3t &v) const
    {
      return vector3t(coo[0]-v.coo[0], coo[1]-v.coo[1], coo[2]-v.coo[2]);
    }

  // 000310
  inline const vector3t operator-(const T &d) const
    {
      return vector3t(coo[0]-d, coo[1]-d, coo[2]-d);
    }

  inline const vector3t operator-(void) const
    {
      return vector3t(-coo[0], -coo[1], -coo[2]);
    }

  //
  // Why can't a non-member operator like this be 'const'?
  //
  inline const friend vector3t operator*(const T &a, const vector3t &v)
    {
      return vector3t(a*v.coo[0], a*v.coo[1], a*v.coo[2]);
    }

  // 030120:
  inline const friend vector3t operator+(const double &a, const vector3t &v)
    {
      return vector3t(a+v.coo[0], a+v.coo[1], a+v.coo[2]);
    }
  inline const friend vector3t operator+(const float &a, const vector3t &v)
    {
      return vector3t(a+v.coo[0], a+v.coo[1], a+v.coo[2]);
    }

  inline T operator*(const vector3t &v) const
    {
      return (coo[0]*v.coo[0] + coo[1]*v.coo[1] + coo[2]*v.coo[2]);
    }

  inline vector3t operator/(const vector3t &v) const
    {
      return (vector3t(coo[1]*v.coo[2]-coo[2]*v.coo[1],
		      coo[2]*v.coo[0]-coo[0]*v.coo[2],
		      coo[0]*v.coo[1]-coo[1]*v.coo[0]));
    }

  //
  // Why can't I make this a const function?
  //
  inline friend T cosangle(const vector3t &v0, const vector3t &v1)
    {
      return ((v0*v1)/sqrt((v0*v0)*(v1*v1)));
    }

  inline const vector3t &clamp(const vector3t &v0, const vector3t &v1)
    {
      coo[0]=std::max(v0.coo[0], std::min(v1.coo[0], coo[0]));
      coo[1]=std::max(v0.coo[1], std::min(v1.coo[1], coo[1]));
      coo[2]=std::max(v0.coo[2], std::min(v1.coo[2], coo[2]));
      return *this;
    }

  //
  // Added 021010:
  // (Hmm... What is the sensible thing to do in such a case,
  // return a new object, like for 'clamp' above, or to just modify
  // *this, like in this definition??)
  //
#ifdef MICROSOFT
  //
  // 030208: I simply don't know where these get defined as macros...
  //         (But I suspect some Microsoft .h-file does it...)
  //
#  undef min
#  undef max
#endif
  inline void min(const vector3t &v)
    {
      coo[0]=std::min(v.x(), x());
      coo[1]=std::min(v.y(), y());
      coo[2]=std::min(v.z(), z());
    }
  inline void max(const vector3t &v)
    {
      coo[0]=std::max(v.x(), x());
      coo[1]=std::max(v.y(), y());
      coo[2]=std::max(v.z(), z());
    }


  //
  // When I do this, how come that I'm allowed to do something like this:
  //
  // vector3t x, y, z, w;
  // z=x.conv(y);
  // z+=w;
  //
  // Aha. Because it's the z I modify, not the return value of conv.
  // This shouldn't work:
  //
  // (x.conv(y)).conv(w)
  //
  // And it does not. Therefore, 'inline vec...' instead of
  // 'inline const vec...'.
  // 
  inline vector3t &conv(const vector3t &v)
    {
      coo[0]*=v.coo[0];
      coo[1]*=v.coo[1];
      coo[2]*=v.coo[2];
      return *this;
    }
  // 050122: Added this variation.
  inline friend vector3t conv2(const vector3t &v1, const vector3t &v2)
    {
      return vector3t(v1.x()*v2.x(), v1.y()*v2.y(), v1.z()*v2.z());
    }

  // 050122: Added this... Note, no checking on components being zero!
  //         So, to scale V to the box with corners A and B:
  //         conv2((B-A).reciprocal(), V-A);
  inline vector3t reciprocal(void) const
    {
      return vector3t(1.0/coo[0], 1.0/coo[1], 1.0/coo[2]);
    }

  // 050122: Might as well add this, then, for better clarity in calling code.
  //         Hmm... 'transform' or something would have been a better name.
  inline vector3t &rescale(const vector3t &mi, const vector3t &ma,
			   const vector3t &new_mi, const vector3t &new_ma)
    {
      *this=conv2(conv2((ma-mi).reciprocal(), *this-mi), new_ma-new_mi)+new_mi;
      return *this;
    }

  inline T length_squared(void) const
    {
      return (coo[0]*coo[0]+coo[1]*coo[1]+coo[2]*coo[2]);
    }

  inline T length(void) const
    {
      return (sqrt(length_squared()));
      //return (sqrt(coo[0]*coo[0]+coo[1]*coo[1]+coo[2]*coo[2]));
    }

  inline void normalize(void)
    {
      const T tmp=1.0/length();
      coo[0]*=tmp;
      coo[1]*=tmp;
      coo[2]*=tmp;
    }

  inline vector3t normalized(void) const
    {
      return (1.0/length())*(*this);
    }

  void print(void) const
    {
      printf("vector3t(%f %f %f)\n",
	     (double)coo[0], (double)coo[1], (double)coo[2]);
    }

  // 041212: Added version which does not print a newline.
  void print2(void) const
    {
      printf("vector3t(%f %f %f)",
	     (double)coo[0], (double)coo[1], (double)coo[2]);
    }

  // 050214: Added version which does not print a newline, fixed size.
  void print3(void) const
    {
      printf("(%7.3f %7.3f %7.3f)",
	     (double)coo[0], (double)coo[1], (double)coo[2]);
    }

  inline T min_coo(void) const
    {
      return std::min(coo[0], std::min(coo[1], coo[2]));
    }

  inline T max_coo(void) const
    {
      return std::max(coo[0], std::max(coo[1], coo[2]));
    }
  
  inline bool dequal(const vector3t<T> &v) const
    {
      return (DEQUAL(x(), v.x()) && DEQUAL(y(), v.y()) && DEQUAL(z(), v.z()));
    }

  inline bool dequal2(const vector3t<T> &v) const
    {
      return (DEQUAL2(x(), v.x()) &&
	      DEQUAL2(y(), v.y()) &&
	      DEQUAL2(z(), v.z()));
    }

  // 030711: For floats
  inline bool dequal3(const vector3t<T> &v) const
    {
      return (DEQUAL3(x(), v.x()) &&
	      DEQUAL3(y(), v.y()) &&
	      DEQUAL3(z(), v.z()));
    }

  // 030209: Useful because we avoid casting to a vector2t.
  inline bool dequal_2d(const vector3t<T> &v) const
    {
      return (DEQUAL(x(), v.x()) && DEQUAL(y(), v.y()));
    }

  // 030709: Rotation in the xy-plane, i.e., around the z-axis.
  inline void rotate_xy(const double cosa, const double sina)
    {
      double oldx=coo[0], oldy=coo[1];
      coo[0] =  cosa*oldx + sina*oldy;
      coo[1] = -sina*oldx + cosa*oldy;
    }

  // 030709: Rotation in the xz-plane, i.e., around the y-axis, but note the
  //         "opposite" orientation...
  inline void rotate_xz(const double cosa, const double sina)
    {
      double oldx=coo[0], oldz=coo[2];
      coo[0] =  cosa*oldx + sina*oldz;
      coo[2] = -sina*oldx + cosa*oldz;
    }

  // 030819: Rotation in the yz-plane, i.e., around the x-axis.
  inline void rotate_yz(const double cosa, const double sina)
    {
      double oldy=coo[1], oldz=coo[2];
      coo[1] =  cosa*oldy + sina*oldz;
      coo[2] = -sina*oldy + cosa*oldz;
    }
};






#endif
