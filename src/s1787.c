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
 * $Id: s1787.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1787

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1787(SISLSurf *ps,double alevel,double aepsge,double epar[],
	   double gpar1[],double gpar2[],int *jstat)
#else
void s1787(ps,alevel,aepsge,epar,gpar1,gpar2,jstat)
     SISLSurf   *ps;
     double alevel;
     double aepsge;
     double epar[];
     double gpar1[];
     double gpar2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : In surface-point intersection in dimention one.
*              To search for endpoints on an intersection curve.
*              The intersection curve must go through the parameter 
*              pair epar. If the function find a closed curve
*              it will return a status message 2. 
*              Else if epar is an edge point the function
*              will return parameter values for the other end point
*              on the intersection curve in gpar1 and a status message
*              1. If epar is an internal point the function will
*              return the parameter values for the two endpoints on
*              the intersection curve in gpar1 and gpar2 and a status
*              message 0.
*
*
* INPUT      : ps    - The surface in intersection.
*              alevel   - The constant value the surface is intersected with.
*              aepsge   - Geometry resolution.
*              epar[2]  - Parameter values for the  intersection point.
*
*
*
* OUTPUT     : gpar1[2] - Parameter values for the first intersection point.
*                         One of the ends are within computer tolerance from
*                         epar, then this point is put in this variable and
*                         the status 1 is returned.
*              gpar2[2] - Parameter values for the second intersection point
*                         on the edges. If closed curve found with no singular
*                         point this point is equal to the first point.
*              jstat   -  status messages  
*                   < 0   : Error.
*                   = 0   : The marching did not succeed.
*                   = 11  : epar == gpar1. Open curve.
*                   = 12  : epar == gpar1. Open curve. gpar2 inside.
*                   = 13  : epar == gpar1. Open curve. gpar1 inside.
*                   = 14  : epar == gpar1. Open curve. Both ends inside.
*                   = 16  : epar == gpar1. Closed curve. gpar2 == gpar1
*                   = 17  : epar == gpar1. Closed curve. gpar2 singular. 
*                   = 21  : epar != output points. Open curve.
*                   = 22  : epar != output points. Open curve. gpar2 inside.
*                   = 24  : epar != output points. Open curve.
*                                                  Output points inside.
*                   = 26  : epar != output points. Closed curve. gpar2 == gpar1
*                   = 27  : epar != output points. Closed curve.
*                                                  gpar2 singular. 
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. July 1989
*
*********************************************************************
*/
{
  int kdeg=1;           /* Indicate that a plane is used               */
  int kk1,kk2,kn1,kn2;  /* Orders and numbers of vertices              */
  int kstat;            /* Status variable                             */
  int kkm1,kkm2;        /* Orders minus 1                              */
  int kincre;           /* Number of doubles in first vertex direction */
  int kpos=0;           /* Position of error                           */
  int ki,kj,kl,kstop;
  int kcur,kgraph;      /* Indicators telling to control type of output
			   from marching                               */
  int kmark1,kmark2,kclose,kmatch1,kmatch2; /* Flags                   */
  double tmax;          /* Box size                                    */
  double tstart,tlength;/* Variables used in Marsdens identity         */
  double tfak;
  double tdum1;         /* Max knot value used in DEQUAL comparing.    */
  double tdum2;         /* Max knot value used in DEQUAL comparing.    */
  double tsum,*sp,*sq;
  double simpli[4];     /* Description of plane                        */
  double *st1,*st2,*scoef; /* Knots and vertices of input surface      */
  double *s3coef=SISL_NULL;  /* 3-D coeff                                   */
			   
  double tepsco = REL_COMP_RES;
  double tepsge;   
  double sval1[2];      /* Limits of parameter plane in first SISLdir      */
  double sval2[2];      /* Limits of parameter plane in second SISLdir     */
  double *spar1,*spar2; /* Pointers to arrays                          */
  double *spar=SISL_NULL;    /* Pointer to allocated values for parameter values*/
  SISLSurf *qs=SISL_NULL;    /* 3-D version of surface                     */
  SISLCurve *qcrv;          /* Curve in parameter plane                   */
  SISLIntcurve *qintcr=SISL_NULL;/* Intersection curve object            */
  kk1   = ps -> ik1;
  kk2   = ps -> ik2;
  kn1   = ps -> in1;
  kn2   = ps -> in2;
  st1   = ps -> et1;
  st2   = ps -> et2;
  scoef = ps -> ecoef;
  sval1[0] = st1[kk1-1];
  sval1[1] = st1[kn1];
  sval2[0] = st2[kk2-1];
  sval2[1] = st2[kn2];
  
  /* Allocate array for 3-D representation of surface */
  
  if((s3coef = newarray(kn1*kn2*3,DOUBLE)) == SISL_NULL) goto err101;
  
  sh1992su(ps,0,aepsge,&kstat);
  if (kstat < 0) goto error;
  
  tmax = ps->pbox->e2max[0][0] - ps->pbox->e2min[0][0];
  
  /* Make description of plane */
  
  simpli[0] = (double)0.0;
  simpli[1] = (double)0.0;
  simpli[2] = (double)1.0;
  simpli[3] = -alevel;
  
  /* Make 3-D description of the surface */
  
  
  /* Make representation of coefficients from Marsdens identity for the
   * function f(t) = t, with the knot vector in first parameter direction
   * scaled to [0,tmax]. This will be used as the x-coordinate in the 3-D
   * representation */
  
  tstart = st1[kk1-1];
  tlength = st1[kn1] - tstart;
  tfak = tmax/tlength;
  kkm1 = kk1 - 1;
  kincre = 3*kn1;
  
  for (ki=0,kl=0,sp=s3coef ; ki<kn1 ; ki++,kl+=3,sp+=3)
    {
      tsum = (double)0.0;
      kstop = ki+kk1;
      for (kj=ki+1;kj<kstop;kj++)
        tsum +=st1[kj];
      
      tsum = (tsum/kkm1-tstart)*tfak;
      
      
      /* Copy x-coordinate to the other vertex rows */
      /*UJK,changed from kj<kn to kj<kn2.*/      
      for (kj=0,sq=sp ; kj<kn2 ; kj++,sq+=kincre) *sq = tsum;
      
    }
  
  /* Make representation of coefficients from Marsdens identity for the
   * function f(t) = t, with the knot vector in second parameter direction
   * scaled to [0,tfak].  This will be used as the x-coordinate in the 3-D
   * representation */
  
  kkm2 = kk2 - 1;
  tstart = st2[kk2-1];
  tlength = st2[kn2] - tstart;
  tfak = tmax/tlength;
  for (ki=0,sp=s3coef+1 ; ki< kn2 ; ki++)
    {
      tsum = (double)0.0;
      kstop = ki+kk2;
      for (kj=ki+1;kj<kstop;kj++)
        tsum +=st2[kj];
      
      tsum  = (tsum/kkm2-tstart)*tfak;
      
      /*  Copy to remaining y-coordinates in first vertex row */
      
      for (kj=0 ; kj<kn1 ; kj++,sp+=3) *sp = tsum;
      
    }
  
  /* Copy z-coordinates */
  
  for (kj=0,sp=s3coef+2,sq=scoef ; kj < kn2 ; kj++)
    for (ki=0 ; ki<kn1 ; ki++,sp+=3,sq++)
      *sp = *sq;
  
  /* Make 3-D surface */
  
  if((qs = newSurf(kn1,kn2,kk1,kk2,st1,st2,s3coef,1,3,1)) == SISL_NULL) goto err101;
  
  kgraph = 0;
  kcur   = 3;

  /* Make an intersection curve object with the parameter value */
  
  if ((spar=newarray(2,DOUBLE))==SISL_NULL) goto err101;
  memcopy(spar,epar,2,DOUBLE);
  
  if((qintcr = newIntcurve(1,2,0,spar,SISL_NULL,0)) == SISL_NULL) goto err101;
  
  kcur = 2;
  kgraph = 0;
  tepsge = tmax*(double)0.01;
  s1313(qs,simpli,kdeg,tepsco,tepsge,tmax,qintcr,kcur,kgraph,&kstat);
  if (kstat==-185) goto war00;
  if (kstat<0) goto error;
  
  /* Identify first and last parameter pair in the intersection curve */
  
  qcrv = qintcr -> ppar1;
  if (qcrv == SISL_NULL) goto war00;
  
  spar1 = qcrv -> ecoef;
  spar2 = spar1 + 2*(qcrv->in)-2;
  /* Check if any of the points lie on the boundary */
  
  kmark1 = 0;
  tdum1 = (double)2.0*max(fabs(sval1[0]),fabs(sval1[1]));
  tdum2 = (double)2.0*max(fabs(sval2[0]),fabs(sval2[1]));

  if (DEQUAL(spar1[0]+tdum1,sval1[0]+tdum1) || DEQUAL(spar1[0]+tdum1,sval1[1]+tdum1) ||
      DEQUAL(spar1[1]+tdum2,sval2[0]+tdum2) || DEQUAL(spar1[1]+tdum2,sval2[1]+tdum2) )
    kmark1 = 1;
  
  kmark2 = 0;

  if (DEQUAL(spar2[0]+tdum1,sval1[0]+tdum1) || DEQUAL(spar2[0]+tdum1,sval1[1]+tdum1) ||
      DEQUAL(spar2[1]+tdum2,sval2[0]+tdum2) || DEQUAL(spar2[1]+tdum2,sval2[1]+tdum2) )
    kmark2 = 1;
  
  /* Check if closed */
  
  kclose = 0;
  if (spar1[0] == spar2[0] && spar1[1] == spar2[1]) kclose = 1;
  
  /* Check if first points matches start point */
  
  kmatch1 = 0;
  if (DEQUAL(epar[0]+tdum1,spar1[0]+tdum1) && DEQUAL(epar[1]+tdum2,spar1[1]+tdum2) ) 
    kmatch1 = 1;
  
  /* Check if second points matches start point */
  
  kmatch2 = 0;
  if (DEQUAL(epar[0]+tdum1,spar2[0]+tdum1) && DEQUAL(epar[1]+tdum2,spar2[1]+tdum2) ) 
    kmatch2 = 1;
  
  /* Check if any point matches start point */
  
  if (kmatch1 == 1 || kmatch2 == 1)
    {
      /*  Start point matches one of the end points, status values in
	  the range 11-19*/
      
      if (kmark1 == 1 && kmark2 == 1 && kclose == 0)
        {
	  /* Open curve, status 11 */
	  *jstat = 11;
	  if(kmatch1==1)
            goto copy;
	  else
            goto invcopy;
        }
      else if (kmark1 ==1 || (kmark2 == 1 && kclose == 0))
	{
	  /* Open curve one point inside status 12 or 13 */
	  
	  if (kmark1 == 1 && kmatch1 == 1)
	    {
	      *jstat = 12;
	      goto copy;
	    }
	  else if (kmark2 == 1 && kmatch2 == 1)
	    {
	      *jstat = 12;
	      goto invcopy;
	    }
	  if (kmark1 == 1 && kmatch2 == 1)
	    {
	      *jstat = 13;
	      goto invcopy;
	    }
	  if (kmark2 == 1 && kmatch1 == 1)
	    {
	      *jstat = 13;
	      goto copy;
	    }
        }
      else if (kclose == 0)
	{
	  /* Both ends inside */
	  *jstat = 14;
	  if(kmatch1==1)
            goto copy;
	  else
            goto invcopy;
	}
      else if(kmatch1 == 1)
	{
	  /* Closed curve, no singularity */
	  *jstat = 16;
	  memcopy(gpar1,spar1,2,DOUBLE);
	  memcopy(gpar2,spar1,2,DOUBLE);
	  goto out;
	}
      else
	{
	  /* Closed curve, with singularity */
	  *jstat=17;
	  memcopy(gpar1,epar ,2,DOUBLE);
	  memcopy(gpar2,spar1,2,DOUBLE);
	  goto out;
	}
    }
  else
    {
      /* epar does not match produced end points, status messages in
	 21-29 the range  */
      
      if (kmark1 ==1 && kmark2 ==1 && kclose == 0)
        {
	  /* Open curve, status 11 */
	  *jstat = 21;
	  memcopy(gpar1,spar1,2,DOUBLE);
	  memcopy(gpar2,spar2,2,DOUBLE);
	  goto out;
        }
      else if (kmark1 ==1 && kclose == 0)
	{
	  /* Open curve one point inside status 12 */
	  *jstat=22;
	  goto copy;
	}
      else if (kmark2 ==1 && kclose == 0)
	{
	  /* Open curve one point inside status 12 */
	  *jstat=22;
	  goto invcopy;
	}
      else if (kclose == 0)
	{
	  /* Both ends inside */
	  *jstat=24;
	  goto copy;
	}
      else if(kmatch1 == 1)
	{
	  /* Closed curve, no singularity */
	  *jstat=26;
	  memcopy(gpar1,spar1,2,DOUBLE);
	  memcopy(gpar2,spar1,2,DOUBLE);
	}
      else
	{
	  /* Closed curve, with singularity */
	  *jstat = 27;
	  memcopy(gpar1,epar ,2,DOUBLE);
	  memcopy(gpar2,spar1,2,DOUBLE);
	  goto out;
	}
    }
  /* Marching produced no curve */
  
 war00: *jstat = 0;
  memcopy(gpar1,epar,2,DOUBLE);
  memcopy(gpar2,epar,2,DOUBLE);
  goto out;
  
 copy:
  memcopy(gpar1,spar1,2,DOUBLE);
  memcopy(gpar2,spar2,2,DOUBLE);
  goto out;
  
 invcopy:
  memcopy(gpar1,spar2,2,DOUBLE);
  memcopy(gpar2,spar1,2,DOUBLE);
  goto out;
  
  /* Error in space allocation */
 err101: 
  *jstat = -101;
  s6err("s1787",*jstat,kpos);
  goto out;
  
  /* Error in lower level function */
 error:
  *jstat = kstat;
  s6err("s1787",*jstat,kpos);
  goto out;
  
 out:
  if (s3coef != SISL_NULL) freearray(s3coef);
  if (qs     != SISL_NULL) freeSurf (qs);
  if (qintcr != SISL_NULL) freeIntcurve(qintcr);
}
