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


#include "sisl-copyright.h"

/*
 *
 * $$
 *
 */
#define S1518

#include "sislP.h"


#if defined(SISLNEEDPROTOTYPES)
void 
s1518(SISLSurf *surf, double point[], double dir[], double epsge,
	   double start[], double end[], double parin[], double parout[],
	   int *stat)
#else
void s1518(surf, point, dir, epsge, start, end, parin, parout, stat)
     SISLSurf   *surf;
     double point[];
     double dir[];
     double epsge;
     double start[];
     double end[];
     double parin[];
     double parout[];
     int    *stat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Newton iteration on the intersection between
*              a 3D NURBS surface and a line.
*              If a good initial guess is given, the intersection will
*              be found quickly. However if a bad initial guess is given,
*              the iteration might not converge.
*              We only search in the rectangular subdomain specified
*              by "start" and "end". This can be the whole domain if desired.
*
*
* INPUT      : surf    - The NURBS surface.
*              point   - A point on the line.
*              dir     - The vector direction of the line
*                        (not necessarily normalized).
*              epsge   - Geometric resolution.
*              start   - Lower limits of search rectangle (umin, vmin).
*              end     - Upper limits of search rectangle (umax, vmax).
*              parin   - Initial guess (u0,v0) for parameter point of
*                        intersection (which should be inside the
*                        search rectangle).
*
*
*
* OUTPUT     : parout  - Parameter point (u,v) of intersection.
*              jstat   - status messages  
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Newton iteration: we assume the surface is close to
*              the tangent plane through the current estimate s(u_0,v_0)
*              and use the intersection of the line with that (parametric)
*              plane to make the next estimate (u_1,v_1).
*
*
* REFERENCES :
*
*
* WRITTEN BY : Michael Floater, SINTEF Oslo, September 2001.
*
*********************************************************************
*/                       
{ 
  int i,j;                  /* For loops. */
  int kstat = 0;            /* Local status variable. */
  int kpos = 0;             /* Error indicator. */
  int num_its = 30;         /* Max number of Newton iterations allowed. */
  double norm1[3];          /* 1st vector normal to dir. */
  double norm2[3];          /* 2nd vector normal to dir. */
  int kleft1=0;             /* Index of knot interval. */
  int kleft2=0;             /* Index of knot interval. */
  int nder=1;               /* evaluate up to derivs of order 1. */
  double parpoint[3];       /* current parameter point. */
  double der[9];            /* surface position and two partial derivatives. */
  double norm[3];           /* surface normal. */
  double a11,a12,a21,a22;   /* Elements of matrix A. */
  double b1,b2;             /* Elements of column vector b. */
  double det;               /* Determinant of A. */
  double du,dv;             /* Increment in parameter point. */
  double su[3],sv[3];       /* two partial derivatives of surface. */
  double pminuss[3];        /* p - s. */
  double epsge2;            /* square of epsge. */

  if(surf->idim != 3) { kstat = -1; goto error; }
  if(start[0] < surf->et1[surf->ik1 - 1]) { kstat = -1; goto error; }
  if(end[0] > surf->et1[surf->in1]) { kstat = -1; goto error; }
  if(start[1] < surf->et2[surf->ik2 - 1]) { kstat = -1; goto error; }
  if(end[1] > surf->et2[surf->in2]) { kstat = -1; goto error; }
  if(parin[0] < start[0] || parin[0] > end[0]) { kstat = -1; goto error; }
  if(parin[1] < start[1] || parin[1] > end[1]) { kstat = -1; goto error; }

  /* Represent line as intersection of two planes, i.e. find
     two vectors norm1 and norm2 of length one,
     perpendicular to the line. */

  s6twonorm(dir,norm1,norm2,&kstat);
    if(kstat < 0) goto error;

  /* printf("norm1 = %lf %lf %lf\n",norm1[0],norm1[1],norm1[2]); */
  /* printf("norm2 = %lf %lf %lf\n",norm2[0],norm2[1],norm2[2]); */

  epsge2 = epsge * epsge;

  parpoint[0] = parin[0];
  parpoint[1] = parin[1];

  /* printf("parpoint = %lf %lf\n",parpoint[0],parpoint[1]); */

  for(i=0; i< num_its; i++)
  {
    /* Evaluate position and 1st derivatives of surface */
  
    s1421(surf,nder,parpoint,&kleft1,&kleft2,der,norm,&kstat);
    if (kstat < 0) goto error;

    /* printf("pos = %lf %lf %lf\n",der[0],der[1],der[2]); */
    /* printf("s_u = %lf %lf %lf\n",der[3],der[4],der[5]); */
    /* printf("s_v = %lf %lf %lf\n",der[6],der[7],der[8]); */
    /* printf("norm = %lf %lf %lf\n",norm[0],norm[1],norm[2]); */
    
    /* We assume that s(u,v) is close to the
       parametric plane s(u_0,v_0) + (u-u_0) * s_u + (v-v_0) * s_v,
       i.e. the tangent plane to s at (u_0,v_0).
       We then find (u_1,v_1) where this plane intersects the
       given straight line. We have represented the line as
       those points x in 3D such that
      
         (x - p) . n1 = 0   (1)
         (x - p) . n2 = 0   (2)
       
       where p is "point" and n1 = "norm1", "n2 = norm2".
       Thus we solve
       
         (s(u_0,v_0) + (u_1-u_0) * s_u + (v_1-v_0) * s_v - p) . n1 = 0   (3)
         (s(u_0,v_0) + (u_1-u_0) * s_u + (v_1-v_0) * s_v - p) . n2 = 0   (4)
    
       which is a 2*2 linear system in the unknowns u_1,v_1.
       This can be written as
  
         A u = b
  
       where A = | s_u . n1   s_v . n1 |    b =  | (p - s) . n1 |
                 | s_u . n2   s_v . n2 |         | (p - s) . n2 |
  
       and u = | du |
               | dv |
    
       and where du = u_1 - u_0, dv = v_1 - v_0.
       We solve the system for du, dv and find u_1 and v_1 afterwards.

       Note that the distance of the point s(u0,v0) to
       the line is precisely the Euclidean norm of the
       right hand side b. Thus if this is within the
       geometric tolerance, we can stop iterating.
    */

    for(j=0; j<3; j++)
    {
      su[j] = der[j+3];
      sv[j] = der[j+6];
      pminuss[j] = point[j] - der[j];
    }
  
    b1 = s6scpr(pminuss,norm1,3);
    b2 = s6scpr(pminuss,norm2,3);

    if(b1 * b1 + b2 * b2 <= epsge2) break;

    a11 = s6scpr(su,norm1,3);
    a12 = s6scpr(sv,norm1,3);
    a21 = s6scpr(su,norm2,3);
    a22 = s6scpr(sv,norm2,3);
  
    det = a11 * a22 - a21 * a12;
    du =  (b1 * a22 - b2 * a12) / det;
    dv =  (a11 * b2 - a21 * b1) / det;
  
    /* Having now found the increments du,dv, update
       the current parameter point. */
  
    parpoint[0] += du;    /* u1 = u0 + du; */
    parpoint[1] += dv;    /* v1 = v0 + dv; */

    /* printf("parpoint = %lf %lf\n",parpoint[0],parpoint[1]); */

   if(surf->cuopen_1 == 1)
   {
     if(parpoint[0] < start[0]) parpoint[0] = start[0];
     if(parpoint[0] > end[0]) parpoint[0] = end[0];
   }
   else // closed in u direction
   {
     if(parpoint[0] < start[0]) parpoint[0] = end[0];
     if(parpoint[0] > end[0]) parpoint[0] = start[0];
   }

   if(surf->cuopen_2 == 1)
   {
   if(parpoint[1] < start[1]) parpoint[1] = start[1];
   if(parpoint[1] > end[1]) parpoint[1] = end[1];
   }
   else // closed in v direction
   {
     if(parpoint[1] < start[1]) parpoint[1] = end[1];
     if(parpoint[1] > end[1]) parpoint[1] = start[1];
   }

  }

  *stat = 1;
  parout[0] = parpoint[0];
  parout[1] = parpoint[1];

  /* printf("parout = %lf %lf\n\n",parout[0],parout[1]); */

  goto out;

  /* Error in lower level routine.  */

  error :
  *stat = kstat;
  s6err("s1518", *stat, kpos);
  goto out;

  out:    

  return;
}

