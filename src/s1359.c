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
 * $Id: s1359.c,v 1.3 2001-03-19 15:58:47 afr Exp $
 *
 */
#define S1359

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1359(double egeo[],double aepsge,int idim,int inbinf,
	   int ipar,double epar[],SISLCurve **rcurve,int *jstat)
#else
void s1359(egeo,aepsge,idim,inbinf,ipar,epar,rcurve,jstat)
     double egeo[];
     double aepsge;
     int    idim;
     int    inbinf;
     int    ipar;
     double epar[];
     SISLCurve  **rcurve;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To represent the curve described in egeo as
*              an Hermit curve on a B-spline format.It hte first and last points
*              have exactly the same position a cyclic basis is made by moving
*              the first and last knots.
*
*
* INPUT      : egeo   - The geometry of the point to be interpolated
*                       The sequence of the information for each point
*                       is: position, unit tangent, curvature vector
*                           and radius of curvature.
*                       When the dimension is 2 this is 7 doubles
*                       When the dimension is 3 this is 10 doubles
*                       Total size of egeo is thus:
*                        idim=2 :  7*inbinf doubles
*                        idim=3 : 10*inbinf doubles
*              aepsge - Geoemtry resolution. Here used for minimal knot
*                       distance.
*              idim   - Dimension of the spcae the points lie in
*                       only 2 and 3 is legal
*              inbinf - Number of points
*              ipar   - Array telling if input parametrization (in epar)
*                       is to be used:
*                        ipar = 0 : Don't use input parametrization
*                        ipar = 1 : Use input parametrization
*
* INPUT/OUTPUT:
*              epar   - Parametrization of the points. ipar determines
*                       if this is input or output
*
* OUTPUT     : rcurve - The curve produced
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, OSLO, Norway, 30. June 1988
*              UJK and VSK 24.10.90, Changed test for recognizing zero curvature
*              Tor Dokken, SI, Cyclic behaviour when closed
*********************************************************************
*/
{
  int kn;             /* Number of vertices                                 */
  int kk = 4;         /* Order of b-spline basis                            */
  int knt;            /* Number of knots produced so far                    */
  int kvert;          /* Pointer to first free variable in vertex array     */
  int kpos =1;        /* Position of error                                  */
  int kstat;          /* Local status variable                              */
  int ki,kj;          /* Running variables in loop                          */
  int kv1,kv2,kv3;    /* Running variables in loop                          */
  int kincre;         /* Number of doubles for each point in egeo           */
  int kcycpos;        /* Indicator telling if cyclic or open geometry       */
  double *sprevp;     /* Pointer to position at start of current segment    */
  double *sprevt;     /* Pointer to tangent  at start of current segment    */
  double *sprevc;     /* Pointer to curvature at start of current segment   */
  double *sprevr;     /* Pointer to radius of curvature start current segment*/
  double snprevt[3];  /* Nomralized version of sprevc                        */
  double *scurp;      /* Pointer to position at end   of current segment     */
  double *scurt;      /* Pointer to tangent  at end   of current segment     */
  double *scurc;      /* Pointer to curvature at end   of current segment    */
  double *scurr;      /* Pointer to radius of curvature end   current segment*/
  double sncurt[3];   /* Normalized version of scurc                         */
  double tcos;        /* Description of angle                                */
  double tl1,tl2;     /* Tangent lengths                                     */
  double tangle;      /* Arcus cosinus if tcos                               */
  double tdist;       /* Distance between start and end of current segment   */
  double tpar;        /* Parameter value at end of segment                   */
  double *st = SISL_NULL;  /* Pointer to knot vector                              */
  double *scoef=SISL_NULL; /* Pointer to vertices                                 */
  double tmpval=aepsge;/* Maximal difference in x, y and z coordinate        */
  double max_dist;
  /* Allocate space for knots and vertices */
  
  if (idim != 2 && idim != 3) goto err105;

  /* No. of points given must be > 1 since it makes no sence to interpolate
     over 0 or 1 point, so check inbinf.*/

  if (inbinf < 2) goto err181;
  
  
  if (idim==2)
    kincre = 7;
  else
    kincre = 10;
  
  
  /* To make sure that we don't make a too long jump in parametrization, find
     maximal difference in x, y and z coordinate. */
  
  if (ipar==0)
    {
      double tmin,tmax,*sp;
      
      for (kj=0 ; kj<idim ; kj++)
        {
	  tmin = *(egeo+kj);
	  tmax = tmin;
	  for (ki=0,sp=egeo+kj ; ki < inbinf ; ki++,sp+=kincre)
            {
	      tmin = MIN(tmin,*sp);
	      tmax = MAX(tmax,*sp);
            }
	  tmpval = MAX(tmpval,(tmax-tmin));
        }
    }
  
  kn = 3*(inbinf-1) + 1;
  scoef = newarray(idim*kn,DOUBLE);
  if (scoef == SISL_NULL) goto err101;
  
  st = newarray(kk+kn,DOUBLE);                                                  
  if (st == SISL_NULL) goto err101;
  
  /* Make four first knots */
  if (ipar==0)
    {
      epar[0] = DZERO;
    }
  
  st[0] = epar[0];
  st[1] = epar[0];
  st[2] = epar[0];
  st[3] = epar[0];
  
  /* Make first vertex */
  memcopy(scoef,egeo,idim,DOUBLE);
  
  
  /* Set pointers to start point, tangent, curvature and radius of curvature
   */
  
  sprevp = egeo;
  sprevt = sprevp + idim;
  sprevc = sprevt + idim;
  sprevr = sprevc + idim;
  
  /* Normalize curvature vector at start */
  
  (void)s6norm(sprevt,idim,snprevt,&kstat);
  
  for (ki=1,knt=4,kvert=idim;ki<inbinf;ki++)
    {
      
      /* For each pair of adjacent points in egeo make an Hermit segment */
      
      /* Set pointers position, tangent, curvature and radius of end of
	 current segment segment */
      
      scurp = sprevp + kincre;
      scurt = sprevt + kincre;
      scurc = sprevc + kincre;
      scurr = sprevr + kincre;
      
      /* Normalize curvature vector at end of segment */
      
      (void)s6norm(scurt,idim,sncurt,&kstat);
      
      /* Make cosine of angle between tangent vectors by making the scalar
	 product of the normalized versions of the two vectors */
      
      tcos = s6scpr(snprevt,sncurt,idim);
      
      /* Find the actual angle by making the arcus tangens of this value */
      
      if (tcos >= DZERO)             
	tcos = MIN((double)1.0,tcos);
      else
	tcos = MAX((double)-1.0,tcos);
      
      tangle = fabs(acos(tcos));
      
      if (tangle < ANGULAR_TOLERANCE) tangle = DZERO;
      
      tdist = s6dist(sprevp,scurp,idim);
      
      /* Make tangent length of start of segment */
      
      
/* UJK and VSK 24.10.90 */
/*      if (DEQUAL(tangle,DZERO) || *sprevr < (double)-1.0) */
      if (DEQUAL(tangle,DZERO) || *sprevr < DZERO)
        {
	  /* Parallel tangents or infinit radius of curvature use 1/3 of
	     the distance between the points as tangent length    */
	  tl1 = tdist/(double)3.0;
        }
      else
        {
	  /* Base tangent length on radius of curvature and opening angle */
	  tl1 = s1325(*sprevr,tangle);
        }
      
      /*  Make tangent length of end of segment */
      
      if (DEQUAL(tangle,DZERO) || *scurr < DZERO)
        {
	  /* Parallel tangents or infinit radius of curvature use 1/3 of
	     the distance between the points as tangent length         */
	  tl2 = tdist/(double)3.0;
        }
      else
        {
	  /* Base tangent length on radius of curvature and opening angle */
	  tl2 = s1325(*scurr,tangle);
        }
      
      /* Make sure that the tangent does not explode due to numeric errors,
	 and make a controlled tangent when the radius is zero or almost zero*/
      
      /* ALA 28.10.93  An improved control was nessesary*/

      if (tangle < 0.1)		max_dist = (double)0.35*tdist;
      else if (tangle < 0.35)	max_dist = (double)0.40*tdist;
      else if (tangle < 0.75)	max_dist = (double)0.50*tdist;
      else 			max_dist = (double)0.70*tdist;

      if ( tl1 > max_dist) tl1 = max_dist;
      if ( tl2 > max_dist) tl2 = max_dist;
      
      /* We want to have a parametrization that is as close as possible to an
         arc length parametrization */                                             
      
      
      if (ipar==0)
        {
	  /* Make parametrization of segment by making an average of arc of a
	     circle with radius sprevr and scurr spanning an angle tangle.
	     If one or both radius infinit use the distance between the points
	     */
	  
	  if (DNEQUAL(*sprevr,(double)-1.0) && 
	      DNEQUAL(*scurr,(double)-1.0))
            {
	      tpar = (double)0.5*tangle*(*sprevr+*scurr);
            }
	  else if (DNEQUAL(*sprevr,(double)-1.0) && 
		   DEQUAL(*scurr,(double)-1.0))
            {
	      tpar = (double)0.5*(*sprevr*tangle + tdist);
            }
	  else if (DEQUAL(*sprevr,(double)-1.0) && 
		   DNEQUAL(*scurr,(double)-1.0))
            {
	      tpar = (double)0.5*(tdist + tangle*(*scurr));
            }
	  else
            {
	      tpar =  tdist;
            }
	  
	  tpar = MAX(tpar,tdist);
	  tpar = MAX(tpar,aepsge);
	  
	  /* Make sure that we don't make a parameter inteval greater than
	     the maximal length of a SISLbox around the input points */
	  
          /* tpar = MIN(tpar,tmpval);
             BOH:220793: Start change */
	  
          if (tangle <= PIHALF)
	     tpar = MIN(tpar, ((double)1.1*tdist));
          else
	     tpar = MIN(tpar,tmpval);
	  
          /* BOH: End change. */
	  
	  if (DEQUAL((epar[ki-1]+tpar),epar[ki-1]))
            {
	      tpar = fabs(epar[ki-1])*(double)0.1;
            }
	  
	  if (DEQUAL(tpar,DZERO))
            {
	      tpar = (double)1.0;
            }
	  
	  epar[ki] = epar[ki-1] + tpar;
	  
        }
      
      /*  Make 3 new knots */
      st[knt]   = epar[ki];
      st[knt+1] = epar[ki];
      st[knt+2] = epar[ki];
      
      /*  Make 3 new vertices of segment */
      
      for (kj=0,kv1=kvert,kv2=kv1+idim,kv3=kv2+idim ; kj<idim ;
	   kj++,kv1++,kv2++,kv3++)
        {
	  scoef[kv1] = sprevp[kj] + tl1*sprevt[kj];
	  scoef[kv2] = scurp[kj]  - tl2*scurt[kj];
	  scoef[kv3] = scurp[kj];
        }
      
      /*  Update pointers */
      sprevp = scurp;
      sprevt = scurt;
      sprevc = scurc;
      sprevr = scurr;
      for (kj=0;kj<idim;kj++) snprevt[kj] = sncurt[kj];
      
      /*  Only update number of vertices if epar[ki-1] != epar[ki] */ 
      
      if (DNEQUAL(epar[ki-1],epar[ki]))
        {
	  kvert+=(3*idim);
	  knt+=3;
        }
    }
  
  /* VSK, 07.94. Moved this statement before inserting the last knot.
     Update number of vertices */
  
  kn = kvert/idim;
  
  /* Insert last knot */
  
  st[kn+kk-1] = st[kn+kk-2];
  
  /* Test if cyclic curve */

  for (ki=0, kcycpos=1 ; ki<idim ; ki++)
    if (egeo[ki] != egeo[(inbinf-1)*kincre+ki]) kcycpos=0;
  
  if (kcycpos ==1)
    {
      st[0] -= (st[kn] - st[kn-1]);
      st[kn+kk-1] += (st[kk] - st[kk-1]);
    }
  
  
  /* Make the curve */
  
  kpos = 1;
  *rcurve = SISL_NULL;
  *rcurve = newCurve(kn,kk,st,scoef,1,idim,1);
  if (*rcurve == SISL_NULL) goto err101;
  
  /* Periodicity flag */
  if (kcycpos)
    {
       test_cyclic_knots(st,kn,kk,&kstat);
       if (kstat<0) goto error;
      if (kstat == 2) (*rcurve)->cuopen = SISL_CRV_PERIODIC;
    }
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1359",*jstat,kpos);
  goto out;

  /* Error in lower level.  */
  
 error: *jstat = kstat;
  s6err("s1359",*jstat,kpos);
  goto out;
  
  /* Error in input, negative relative tolerance given */
  
 err105: *jstat = -105;
  s6err("s1359",*jstat,kpos);
  goto out;

  /* Error in input, to few points given(inbinf < 2. */

  err181: *jstat = -181;
   s6err("s1359",*jstat,kpos);
   goto out;
  
  /* Free allocated arrays */
 out:
  
  
  if (st != SISL_NULL)    freearray(st);
  if (scoef != SISL_NULL) freearray(scoef);
  
  
  return;
}
