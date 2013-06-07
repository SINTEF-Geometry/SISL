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
 * $Id: s1788.c,v 1.3 2005-02-28 09:04:49 afr Exp $
 *
 */


#define S1788

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1788(SISLSurf *ps1,SISLSurf *ps2,double aepsge,double epar[],
	   double gpar1[],double gpar2[],int *jstat)
#else
void s1788(ps1,ps2,aepsge,epar,gpar1,gpar2,jstat)
     SISLSurf   *ps1;
     SISLSurf   *ps2;
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
* PURPOSE    : In surface-surface intersection in dimention three.
*              To search for endpoints on an intersection curve.
*              The intersection curve must go through the parameter 
*              values epar. If the function find a closed curve
*              it will return a status message 2. 
*              Else if epar is an edges point the function
*              will return parameter values for the other end point
*              on the intersection curve in gpar1 and a status message
*              1. If epar is an intrnal point the function will
*              return the parameter values for the two endpoints on
*              the intersection curve in gpar1 and gpar2 and a status
*              message 0.
*
*
* INPUT      : ps1      - The first surface in intersection.
*              ps2      - The second surface in intersection.
*              aepsge   - Geometry resolution.
*              epar[4]  - Parameter values for the  intersection point.
*
*
*
* OUTPUT     : gpar1[4] - Parameter values for the first intersection point.
*                         One of the ends are within computer tolerance from
*                         epar, then this point is put in this variable and
*                         the status 1 is returned.
*              gpar2[4] - Parameter values for the second intersection point
*                         on the edges. If closed curve found with no singular
*                         point this point is equal to the first point.
*              jstat   - status messages    
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

*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. July 1989
* REVISED BY : UJK, November 1990 : Changed call to newIntcurve        
*********************************************************************
*/
{
  int kpos=0;           /* Position of error                           */
  int kk1,kk2,kn1,kn2;  /* Orders and numbers of vertices              */
  int kstat;            /* Status variable                             */
  int kmark1,kmark2,kclose,kmatch1,kmatch2; /* Flags                   */
  int kcur,kgraph;      /* Indicators telling to control type of output
			   from marching                               */
  double sval1[2];      /* Limits of parameter plane in first SISLdir      */
  double sval2[2];      /* Limits of parameter plane in second SISLdir     */
  double sval3[2];      /* Limits of parameter plane in third SISLdir      */
  double sval4[2];      /* Limits of parameter plane in fourth SISLdir     */
  double *st1,*st2;     /* Knots and vertices of input surface         */
  double tmax;          /* Box size                                    */
  double tepsge;   
  double *spar11,*spar12;/* Pointers to arrays                         */
  double *spar21,*spar22;/* Pointers to arrays                         */

  /* UJK, Nov 1990 */
  /*  double *spar=SISL_NULL; */
  double *spar1=SISL_NULL;     /* Pointer to allocated values for parameter values*/
  double *spar2=SISL_NULL;     /* Pointer to allocated values for parameter values*/

  SISLCurve *qcrv;           /* Curve in parameter plane                    */
  SISLIntcurve *qintcr=SISL_NULL; /* Intersection curve object            */
  
  /* Find limits of parameter plane */
  
  kk1   = ps1 -> ik1;
  kk2   = ps1 -> ik2;
  kn1   = ps1 -> in1;
  kn2   = ps1 -> in2;
  st1   = ps1 -> et1;
  st2   = ps1 -> et2;
  sval1[0] = st1[kk1-1];
  sval1[1] = st1[kn1];
  sval2[0] = st2[kk2-1];
  sval2[1] = st2[kn2];
  
  kk1   = ps2 -> ik1;
  kk2   = ps2 -> ik2;
  kn1   = ps2 -> in1;
  kn2   = ps2 -> in2;
  st1   = ps2 -> et1;
  st2   = ps2 -> et2;
  sval3[0] = st1[kk1-1];
  sval3[1] = st1[kn1];
  sval4[0] = st2[kk2-1];
  sval4[1] = st2[kn2];
  
  
  /* Make maximal step length based on box-size of surface */
  
  sh1992su(ps1,0,aepsge,&kstat);
  if (kstat < 0) goto error;
  
  tmax = MAX(ps1->pbox->e2max[0][0] - ps1->pbox->e2min[0][0],
	     ps1->pbox->e2max[0][1] - ps1->pbox->e2min[0][1]);
  tmax = MAX(tmax,ps1->pbox->e2max[0][2] - ps1->pbox->e2min[0][2]);
  
  sh1992su(ps2,0,aepsge,&kstat);
  if (kstat < 0) goto error;
  
  tmax = MAX(tmax,ps2->pbox->e2max[0][0] - ps2->pbox->e2min[0][0]);
  tmax = MAX(tmax,ps2->pbox->e2max[0][1] - ps2->pbox->e2min[0][1]);
  tmax = MAX(tmax,ps2->pbox->e2max[0][2] - ps2->pbox->e2min[0][2]);
  
  tepsge = tmax * (double)0.01;
  

  kgraph = 0;
  kcur   = 3;
  
  /* Make an intersection curve object with the parameter value */
  /* UJK, Nov 1990 */  
  /*  if ((spar=newarray(4,DOUBLE))==SISL_NULL) goto err101;
      memcopy(spar,epar,4,DOUBLE); */
  if ((spar1=newarray(2,DOUBLE))==SISL_NULL) goto err101;
  memcopy(spar1,epar,2,DOUBLE);
  if ((spar2=newarray(2,DOUBLE))==SISL_NULL) goto err101;
  memcopy(spar2,epar+2,2,DOUBLE);
  
  /* UJK, Nov 1990 */  
  /* if((qintcr = newIntcurve(1,2,2,epar,epar+2,0)) == SISL_NULL) goto err101; */
  if((qintcr = newIntcurve(1,2,2,spar1,spar2,0)) == SISL_NULL) goto err101;
  
  kcur = 2;
  kgraph = 0;
  tmax = (double)0.0;

  
  s1310(ps1,ps2,qintcr,tepsge,tmax,kcur,kgraph,&kstat);
  
  if (kstat==-185) goto war00;
  if (kstat<0) goto error;
  
  /* Identify first and last parameter pair in the intersection curve */
  
  qcrv = qintcr -> ppar1;
  if (qcrv == SISL_NULL) goto war00;
  
  spar11 = qcrv -> ecoef;
  spar21 = spar11 + 2*(qcrv->in)-2;
  
  
  qcrv = qintcr -> ppar2;
  if (qcrv == SISL_NULL) goto war00;
  
  spar12 = qcrv -> ecoef;
  spar22 = spar12 + 2*(qcrv->in)-2;
  
  /* Check if any of the points lie on the boundary */
  
  kmark1 = 0;
  if (spar11[0] == sval1[0] || spar11[0] == sval1[1] ||
      spar11[1] == sval2[0] || spar11[1] == sval2[1] ||
      spar12[0] == sval3[0] || spar12[0] == sval3[1] ||
      spar12[1] == sval4[0] || spar12[1] == sval4[1] ) kmark1 = 1;
  
  kmark2 = 0;
  if (spar21[0] == sval1[0] || spar21[0] == sval1[1] ||
      spar21[1] == sval2[0] || spar21[1] == sval2[1] ||
      spar22[0] == sval3[0] || spar22[0] == sval3[1] ||
      spar22[1] == sval4[0] || spar22[1] == sval4[1] ) kmark2 = 1;
  
  
  /* Check if closed */
  
  kclose = 0;
  if (spar11[0] == spar21[0] && spar11[1] == spar21[1]  &&
      spar12[0] == spar22[0] && spar12[1] == spar22[1]     ) kclose = 1;
  
  /* Check if first points matches start point */
  
  kmatch1 = 0;
  if (DEQUAL(epar[0],spar11[0]) && DEQUAL(epar[1],spar11[1]) &&
      DEQUAL(epar[2],spar12[0]) && DEQUAL(epar[3],spar12[1])   ) kmatch1 = 1;
  
  /* Check if second points matches start point */
  
  kmatch2 = 0;
  if (DEQUAL(epar[0],spar21[0]) && DEQUAL(epar[1],spar21[1]) &&
      DEQUAL(epar[2],spar22[0]) && DEQUAL(epar[3],spar22[1]) ) kmatch2 = 1;
  
  
  /* Check if any point matches start point */
  
  if (kmatch1 == 1 || kmatch2 == 1)
    {
      /* Start point matches one of the end points, status values in
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
	  memcopy(gpar1,  spar11,2,DOUBLE);
	  memcopy(gpar1+2,spar12,2,DOUBLE);
	  memcopy(gpar2,  gpar1, 4,DOUBLE);
	  goto out;
	}
      else
	{
	  /* Closed curve, with singularity */
	  *jstat=17;
	  memcopy(gpar1,  epar,  4,DOUBLE);
	  memcopy(gpar2,  spar11,2,DOUBLE);
	  memcopy(gpar2+2,spar12,2,DOUBLE);
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
	  goto copy;
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
	  memcopy(gpar1,  spar11,2,DOUBLE);
	  memcopy(gpar1+2,spar12,2,DOUBLE);
	  memcopy(gpar2,  gpar1, 4,DOUBLE);
	}
      else
	{
	  /* Closed curve, with singularity */
	  *jstat = 27;
	  memcopy(gpar1,  epar,  4,DOUBLE);
	  memcopy(gpar2  ,spar11,2,DOUBLE);
	  memcopy(gpar2+2,spar12,2,DOUBLE);
	  goto out;
	}
    }

  /* Marching produced no curve */

 war00: 
  *jstat = 0;
  memcopy(gpar1,epar,4,DOUBLE);
  memcopy(gpar2,epar,4,DOUBLE);
  goto out;

 copy:
  memcopy(gpar1,  spar11,2,DOUBLE);
  memcopy(gpar1+2,spar12,2,DOUBLE);
  memcopy(gpar2,  spar21,2,DOUBLE);
  memcopy(gpar2+2,spar22,2,DOUBLE);
  goto out;

 invcopy:
  memcopy(gpar1,  spar21,2,DOUBLE);
  memcopy(gpar1+2,spar22,2,DOUBLE);
  memcopy(gpar2,  spar11,2,DOUBLE);
  memcopy(gpar2+2,spar12,2,DOUBLE);
  goto out;

  /* Error in space allocation */
 err101: 
  *jstat = -101;
  s6err("s1788",*jstat,kpos);
  goto out;

  /* Error in lower level function */
 error:
  *jstat = kstat;
  s6err("s1788",*jstat,kpos);
  goto out;

 out:
  if (qintcr != SISL_NULL) freeIntcurve(qintcr);
}
