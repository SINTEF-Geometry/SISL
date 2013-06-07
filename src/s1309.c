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
 * $Id: s1309.c,v 1.3 2001-03-19 15:58:43 afr Exp $
 *
 */
#define S1309

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
s1309(double epnt[],double edir[],double eimpli[],int ideg,int *jstat)
#else
double s1309(epnt,edir,eimpli,ideg,jstat)
     double epnt[];
     double edir[];
     double eimpli[];
     int    ideg;
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To find the shortest distance between a point and the
*              projection of the point along a direction vector on to
*              an implicit represented surface. The distance will
*              have a sign, positive if the projection is along the
*              direction of the direction vector, else it is negative.
*              For torus surfaces the projection direction is not
*              used: The closest point is calculated and the distance to
*              this returned.
*              If ideg==1003,1004,1005, the the difference between PI/2 and the
*              angle between the two vectors which whose scalar product
*              wants to be as close as possible to zero, is found.
*
* INPUT      : epnt   - The point to be projected
*              edir   - Projection direction. Not used if ideg==1003,1004,1005
*              eimpli - Description of implicit surface
*              ideg   - Degree of implicit surface
*                        ideg=1:    Plane
*                        ideg=2;    Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*
*
* OUTPUT     : s1309  - The distance, positive if along positive direction
*                       of edir, negative if along negative direction of
*                       edir.
*                       If ideg==1003,1004,1005, the the difference between PI/2 and the
*                       angle between the two vectors which whose scalar product
*                       wants to be as close as possible to zero, is found (in radians).
*              jstat  - status messages
*                       = 0      : ok, iteration converged
*                       < 0      : error
*
*
* METHOD     :
*
* USE        : The function can be used for preparing the calculation
*              of the projected point by adding the product of the distance
*              and the normalized versions of edir to epnt.
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 4-July-1988
* Revised by : Mike Floater, SI, 1991-01
*              Improved the routine for parallel silhouettes (ideg=1003) and
*              added perspective and circular silhouettes (ideg=1004,ideg=1005)
* Changed by : Per OEyvind Hvidsten, SINTEF, 11-94.
*              Added initialization of tcurdst in declaration.
*
*********************************************************************
*/
{
  double sdir[3];         /* Normilized direction vector          */
  double tb1,ta11,ta12;   /* Dummy variables                      */
  double tsum,t1,t2,tdum1;/* Dummy variables                      */
  double tcurdst=0.0;     /* The distance                         */
  double sq[4];           /* Array used for temporary results     */
  int ksize;              /* Number of doubles for storage of derivatives
			     and normal vector */
  int ksizem3;            /* ksize - 3                                      */
  int    kstat;           /* Local status variable                */
  int    kdim=3;          /* Dimesnon of 3-D space                */
  int    ki,kj,kl,kp;     /* Control variables in loop            */
  int    kpos=1;          /* Position of error                    */



  /* If ideg=1,2 or 1001 then only derivatives up to second order
     are calculated, then 18 doubles for derivatives and 3 for the
     normal vector are to be used for calculation of points in the
     spline surface. For ideg=1003,1004,1005 we have a silhouette curve and
     derivatives up to the third are to be calculated,
     thus 30 +3 a total of 33 doubles are to be calculated */

  if (ideg==1003 || ideg==1004 || ideg==1005)
    {
      ksize = 33;
    }
  else
    {
      ksize = 21;
    }

  ksizem3 = ksize -3;
  (void)s6norm(edir,kdim,sdir,&kstat);
  if (kstat < 0) goto error;

  if (ideg==1)
    {
      /* Put parametric representation of projection line into
	 implicit equation of plane */
      tb1  = s6scpr(eimpli,epnt,kdim);
      ta11 = s6scpr(eimpli,sdir,kdim);
      if( ta11 == (double)0.0) goto war02;
      tcurdst = -(eimpli[3]+tb1)/ta11;
    }
  else if (ideg==2)
    {

      /*  Find distance from the new point to the implicit surface, by
	  intersecting  the straight line throught the point and with
	  sdir as normalized direction vector, normalized version of 3
	  first coordinates of sq the implicit surface and P0 = (P,1).

	  This problem can be written:

	  T                T    2                    T
	  P0 A P0  + 2 t P0 A sdir  + t  (sdir,0) A (sdir,0)  = 0

	  T
	  We have to calulate calculate  tb1=P0 A P  and qs=P0 A, thus:
	  T    2                    T
	  tb1 + 2 t sq sdir  + t  (sdir,0) A (sdir,0)  = 0

	  Thus the first step is to calculate qs = A (P,1)
	  */
      for (ki=0;ki<4;ki++)
        {
	  tsum = eimpli[12+ki];
	  for (kj=0,kl=ki ; kj<3 ; kj++,kl+=4)
            {
	      tsum +=(eimpli[kl]*epnt[kj]);
            }
	  sq[ki] = tsum;
        }

      tb1  = s6scpr(epnt,sq,kdim) + sq[3];

      ta11 = (double)2.0*s6scpr(sq,sdir,kdim);

      ta12 = (double)0.0;
      for (ki=0,kl=0;ki<3;ki++,kl+=4)
        {
	  kp = kl;
	  for (kj=0;kj<3;kj++,kp++)
            {
	      ta12 += sdir[ki]*eimpli[kp]*sdir[kj];
            }
        }

      /*  Now our equation system is:
	  2
	  ta12 t  + ta11 t + tb1 = 0

	  we want the root with the smallest absolute value:
	  2
	  t = (-ta11 +/- sqrt(ta11 -4ta12 tb1))/(2ta12)
	  */
      if (DNEQUAL(ta12,(double)0.0))
        {
	  tdum1 = ta11*ta11 - (double)4.0*ta12*tb1;
	  if (tdum1 < DZERO) goto war02;
	  tdum1 = sqrt(tdum1);
	  t1 = (-ta11 + tdum1)/((double)2.0*ta12);
	  t2 = (-ta11 - tdum1)/((double)2.0*ta12);
	  if (fabs(t1)<fabs(t2))
            {
	      tcurdst = t1;
            }
	  else
            {
	      tcurdst = t2;
            }
        }
      else if(DNEQUAL(ta11,(double)0.0))
        {
	  tcurdst = tb1/ta11;
        }
      else
        {
	  /*      Unsolvable system */
	  goto war02;
        }
    }
  else if (ideg==1001)
    {
      /*  Torus surface */

      double *scentr;  /* The center of the torus */
      double *snorm;   /* The normal of the torus symmetry plane */
      double tbigr;    /* The big radius of the torus */
      double tsmalr;   /* The small radius of the torus */
      double sdum1[3]; /* Temporary storage for point */
      double sdum2[3]; /* Temporary storage for point */
      double tproj;    /* Projection of vector onto snorm */


      scentr = eimpli;
      snorm  = eimpli+3;
      tbigr  = *(eimpli+6);
      tsmalr = *(eimpli+7);

      /*  Find projection of vector from torus center on to torus axis */
      s6diff(epnt,scentr,kdim,sdum1);
      tproj = s6scpr(sdum1,snorm,kdim);

      /*  Project vector from torus center to epnt onto torus plane */
      for (ki=0;ki<kdim;ki++)
        sdum2[ki] = sdum1[ki] - tproj*snorm[ki];
      (void)s6norm(sdum2,kdim,sdum2,&kstat);
      if (kstat<0) goto error;

      /*  Find vector from torus circle to epnt */
      for (ki=0;ki<kdim;ki++)
        sdum1[ki] = sdum1[ki] - tbigr*sdum2[ki];

      /*  Find length of this vector and compare with tsmalr */

      tcurdst = fabs(s6length(sdum1,kdim,&kstat)-tsmalr);
      if (kstat<0) goto error;
    }
  else if (ideg==1003)
    {
      /*  Silhouette line/curve */

      double sdum1[3]; /* Temporary storage for point */

      (void)s6norm(epnt+ksizem3,kdim,sdum1,&kstat);
      if (kstat<0) goto error;

      /*  eimpli[0,1,2] is assumed to, be normalized */

      t1 = s6scpr(sdum1,eimpli,kdim);

      /*  t1 now contains the cosine of the angle between the view direction
	  and the normal vector. This is Equal to sin of PI/2 minus the angle
	  between the vectors, thus the actual angle can be calulated */

      tcurdst = asin(t1);
      tcurdst = fabs(tcurdst);
    }

  else if (ideg==1004)
    {
      /*  Perspective silhouette line/curve */

      double sdum1[3],sdum2[3]; /* Temporary storage for point */

      s6diff(epnt,eimpli,kdim,sdum1);
      (void)s6norm(sdum1,kdim,sdum2,&kstat);
      /* OK if sdum1 is zero -- tcurdst will be zero as well  */
      (void)s6norm(epnt+ksizem3,kdim,sdum1,&kstat);


      t1 = s6scpr(sdum1,sdum2,kdim);

      /*  t1 now contains the cosine of the angle between the direction of the
          point epnt relative to the eyepoint E (in eimpli)
	  and the normal vector. This is Equal to sin of PI/2 minus the angle
	  between the vectors, thus the actual angle can be calulated */

      tcurdst = asin(t1);
      tcurdst = fabs(tcurdst);
    }

  else if (ideg==1005)
    {
      /*  Circular silhouette line/curve */

      double sdum1[3],sdum2[3]; /* Temporary storage for point */
      double *bvec=eimpli+3;

      s6diff(epnt,eimpli,kdim,sdum1);
      s6crss(epnt+ksizem3,sdum1,sdum2);
      (void)s6norm(sdum2,kdim,sdum1,&kstat);
      /* OK if sdum2 is zero -- tcurdst will be zero as well  */

      /*  bvec  = eimpli[3,4,5] is assumed to, be normalized */

      t1 = s6scpr(sdum1,bvec,kdim);

      /*  t1 now contains the cosine of the angle between
          (the cross product of the normal vector and the direction of the
          point epnt relative to the point Q (in eimpli))
	  and the direction vector B. This is Equal to sin of PI/2 minus the angle
	  between the vectors, thus the actual angle can be calulated */

      tcurdst = asin(t1);
      tcurdst = fabs(tcurdst);
    }

  *jstat = 0;
  goto out;

  /* Projection not possible */
 war02: *jstat = 2;
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  s6err("s1309",*jstat,kpos);
  goto out;

 out:
  return(tcurdst);
}
