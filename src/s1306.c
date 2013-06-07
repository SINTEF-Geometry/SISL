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
 * $Id: s1306.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */


#define S1306

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1306(double ep[],double eparp[],double eimpli[],int ideg,
	   double egeo3d[],double egeop[],int *jstat)
#else
void s1306(ep,eparp,eimpli,ideg,egeo3d,egeop,jstat)
     double ep[];
     double eparp[];
     double eimpli[];
     int    ideg;
     double egeo3d[];
     double egeop[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the radius of curvature at the intersection
*              point between a B-spline surface and an implicit represented
*              surface. The first and second derivatives
*              of the intersection point in the B-spline surface are give as
*              input. All coordinates are assumed to be in 3-D.
*
*
*
* INPUT      : ep     - 0-2 order derivatives of B-spline surface.
*                       For ideg=1,2 and 1001 the sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) derivative
*                       and normal. (21 numbers)
*                       For ideg=1003,1004,1005 the second derivatives are followed
*                       by the third derivatives and the normal (33 numbers)
*                       Compatible with output of s1421
*              eparp  - Parameter pair in B-spline surface of point
*              eimpli - Description of the implicit surface
*              ideg   - The degree of the implicit surface
*                        ideg=1:    Plane              
*                        ideg=2;    Quadric surface
*                        ideg=1001: Torus surface
*                        ideg=1003: Silhouette line parallel projection
*                        ideg=1004: Silhouette line perspective projection
*                        ideg=1005: Silhouette line circular projection
*
*
* OUTPUT     : 
*              jstat  - status messages
*                         = 11     : Fuzzy singular intersection point found
*                         = 10     : Singular intersection point found
*                         = 2      : Singular intersection point found
*                         = 1      : Curvature radius infinit
*                         = 0      : ok, curvature radius
*                         < 0      : error
*              egeo3d - 3-D geometry description of the intersection. The
*                       array contains: position, unit tangent, curvature
*                       and radius of curvature. (A total of 10 numbers)
*                       A radius of curvature =-1, indicates that the radius
*                       of curvature is infinit.
*              egeop  - Description of the intersection in the parameter plane
*                       of the B-spline surface. The array contains: position,
*                       unit tangent, curvature and radius of curvature.
*                       (A total of 7 numbers)
*
* METHOD     : We put the parametric surface into the implicit equation
*              and get a function f(s,t)=0. By assuming that u and v are
*              function s of a variable u we get:
*
*                f(s(u),y(u)) = 0
*
*              By making the derivative with respect to u we get:
*
*                f s  + f t  = 0
*                 s u    t u    
*
*              By assuming that u=s or u=t, we get one derivative equal
*              to 1, and the other can be calculated. By taking one more
*              derivative with respect to u we get:
*
*                    2                           2
*                f  s  + 2 f  s t  + f s   + f  t  + f t   = 0.
*                 ss u      st u u    s uu    tt u    t uu
*
*              Since s=u or t=u we know that one of the double derivatives
*              is zero and the other can be caluclated.
*
*              Based on these derivatives in the parameter plane the
*              derivatives of the intersection curve with respect to u can
*              be calculated, and futher on the actual curvature of the
*              intersection curve at the input point.
*
*              In the case that both f = 0 and f = 0 the second equation
*                                     s         t
*              can be used for calculating s  and t
*                                           u      u.
*
* USE:         This function is only working in 3-D
*
* REFERENCES : 
*-
* CALLS      : s6scpr,fabs,s1307,s6err
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway, 30 May 1988
* Revised by : Tor Dokken, SI, Oslo, March 1989.
*              Initiating tangent etc. to 0 when singular point found
* Revised by : Mike Floater, SI, 1991-01
*                   Add perspective and circular silhouettes (ideg=1004,ideg=1005)
*
*********************************************************************
*/
{
  int fuzzy_sing = FALSE;
  int sing = FALSE;
  int kstat = 0;          /* Local status variable                       */
  int ki,kj,kl;           /* Control variables in loop                   */
  int kpos = 0;           /* Position of error                           */
  int ksize;              /* Number of doubles for storage of derivateves
			     and normal vector */      
  int ksizem3;            /* ksize - 3                                   */
  double *sps;            /* Pointer to dP/ds                            */
  double *spt;            /* Pointer to dP/dt                            */
  double *spss;           /* Pointer to ddP/(dsds)                       */
  double *spst;           /* Pointer to ddP/(dsdt)                       */
  double *sptt;           /* Pointer to ddP/(dtdt)                       */
  double tfs,tft,tfss;    /* Derivatives of parametric surface put into  */
  double tfst,tftt;       /* the implicit equation                       */
  double tsu,ttu,tsuu,ttuu;/* Derivatives of parameter direction         */
  double sder[6];         /* Derivatives of parametric surface           */
  double snorm[3];        /* Normal vector of implicit surface at ep     */
  double sp[9];           /* Points, first and second deriv. of intcur   */
  
  
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
  
  /* Calculated derivatives of the parametric surface put into the implicit
     surface at the point ep */
  
  s1331(ep,eimpli,ideg,2,sder,snorm,&kstat);
  if (kstat<0) goto error;
  
  tfs  = sder[1];
  tft  = sder[2];
  tfss = sder[3];
  tfst = sder[4];
  tftt = sder[5];
  
  
  /* Calculate ds/du and dt/du */
  /* UJK, aug.92 */
  
  /* if (DEQUAL(tfs,(double)0.0) && 
      DEQUAL(tft,(double)0.0) ) */

  /* UJK, 13.08.93, save suzzy singular information. */
  if (fabs(tfs) < 0.000001 && 
      fabs(tft) < 0.000001) fuzzy_sing = TRUE;
  else
     fuzzy_sing = FALSE;
  
  if (DEQUAL((double)1.0 + tfs,(double)1.0) && 
      DEQUAL((double)1.0 + tft,(double)1.0) )
    {
      double tdum1,tdum2,tafss,tafst,taftt;
      
      /* Singular point found, copy position and derivative value of input to
	 output, if not second order derivatives uniqely describes a tangent*/
      
      memcopy(egeo3d,ep,3,DOUBLE);
      memcopy(egeop,eparp,2,DOUBLE);
      
      tafss = fabs(tfss);
      tafst = fabs(tfst);
      taftt = fabs(tftt);
      
      tdum1 = tfst*tfst - tfss*tftt;
      tdum2 = MAX(tafss,taftt);
      tdum2 = MAX(tafst,tdum2);
      
      for (ki=3 ; ki<10 ; ki++) egeo3d[ki] = DZERO;
      for (ki=2 ; ki<7  ; ki++) egeop[ki]  = DZERO;
      
      if (DEQUAL(tdum2+tdum1,tdum2) &&
	  (DNEQUAL(tafss+tafst,tafst) || DNEQUAL(taftt+tafst,tafst)) )
        {
	  /* A unique tangent can be calculated */
	  
	  if (tafss>taftt)
            {
	      tsu = -tfst/tfss;
	      ttu = 1;
            }
	  else
            {
	      tsu = 1;
	      ttu = -tfst/tftt;
            }
	  tsuu = 0.0;
	  ttuu = 0.0;
	  /* Flag singular case */
	  sing = TRUE;
	}
      else
        {
	  /* A tangent can not be calulated */
	  goto war02;
        }
    }
  else
    {
      /* A noneisngular point found */
      if (fabs(tfs) > fabs(tft))
	{
	  /*  Use u=t */
	  
	  tsu = -tft/tfs;
	  ttu = (double)1.0;   
	  
	  tsuu = -(tfss*tsu*tsu + (double)2.0*tfst*tsu*ttu + tftt*ttu*ttu)/tfs;
	  ttuu = (double)0.0;
	}
      else
	{ 
	  /*  Use u=t */
	  
	  tsu = (double)1.0;
	  ttu = -tfs/tft;  
	  
	  tsuu = (double)0.0;
	  ttuu = -(tfss*tsu*tsu + (double)2.0*tfst*tsu*ttu + tftt*ttu*ttu)/tft;
	}
    }
  
  /* The calculation of the derivatives of one parameter direction with
     respect to the other is dependent on the degree of the implicit equation*/
  
  /* Set local pointers */
  
  sps = ep + 3;
  spt = ep + 6;
  spss = ep + 9;
  spst = ep + 12;
  sptt = ep + 15;
  
  /* Make description of intersection point in 3-D including
     curvature and radius of curvature, first make description of
     intersection point in parameter plane */
  
  /* We will now express the intersection curve locally as a function
   *  of a w-parameter.                                              
   *
   *  c(w) = p(s(w),t(w))
   *
   *  This gives the derivative
   *                         
   *   
   *  c' = P s' + P t'                                                      
   *        s      t
   *
   *  And the second derivative
   *
   *            2                   2
   *  c" = P  s'  + 2P  s't' + P  t' + P s" + P t"
   *        ss        st        tt      s      t  
   */
  
  for (ki=0,kj=3,kl=6 ; ki<3 ; ki++,kj++,kl++)
    {
      
      /*  Copy position */
      
      sp[ki] = ep[ki];
      
      /*  Make tangent  */
      
      sp[kj] = sps[ki]*tsu + spt[ki]*ttu;
      
      /*  Make curvature */
      
      sp[kl] = spss[ki]*tsu*tsu + (double)2.0*spst[ki]*tsu*ttu + 
	sptt[ki]*ttu*ttu + sps[ki]*tsuu + spt[ki]*ttuu;
    }
  
  /* Make 3-D curvature and radius of curvature */
  
  s1307(sp,3,egeo3d,&kstat);
  if (kstat < 0 )goto error;
  
  /* Make description of intersection point in parameter plane including
     curvature and radius of curvature, first make description of
     intersection point in parameter plane */
  
  sp[0] = eparp[0];
  sp[1] = eparp[1];
  sp[2] = tsu;
  sp[3] = ttu;
  sp[4] = tsuu;
  sp[5] = ttuu;
  
  s1307(sp,2,egeop,&kstat);
  if (kstat < 0) goto error;
  
  /* Everyting is ok */
  
  *jstat = 0;
  goto out;
  
  
  /* SISLPoint lying on torus axis */
  
 war02: *jstat=2;
  goto out;
  
  /* Error in lower level function */
 error:
  *jstat = kstat;
  s6err("s1306",*jstat,kpos);
  goto out;
  
 out: 
    if (sing && *jstat>=0) *jstat = 10;
    else if (fuzzy_sing && *jstat>=0) *jstat = 11;
  return;
}
