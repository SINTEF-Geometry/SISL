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
 * $Id: s1612.c,v 1.2 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1612

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1612(SISLCurve *pc,double aepsge,double **gpoint,int *jnbpnt,int *jleng,int *jstat)
#else
void s1612(pc,aepsge,gpoint,jnbpnt,jleng,jstat)
     SISLCurve  *pc;
     double aepsge;
     double **gpoint;
     int    *jnbpnt;
     int    *jleng;
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate a set of points on a B-spline curve. 
*              The straight lines between the points will not deviate
*              more than aepsge from the B-spline curve at any point.
*             
* INPUT      : pc     - The input B-spline curve.   
*              aepsge - Geometry resolution, maximum distance allowed between
*                       the curve and the straight lines to be calculated. 
*
* INPUT/OUTPUT :
*              jnbpnt - No. of calculated points until now 
*              jleng  - no. of allocated doubles in gpoint until now
*              gpoint - Calculated points    
*
* OUTPUT     : jstat  - status messages         
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : This routine recursively calls itself.
*              First the curve is split at kk-1 internal multiple knots.
*              Each curve segment is then treated: 
*              The distance between the control polygon and a straight
*              line from the first to the last vertex is found. If the 
*              distance is grater than aepsge, the curve is split again
*              as follows: 
*              If any internal knots exist the curve is split in the
*              middlemost knot; if not, the curve is split in the middle of
*              the parameter interval. 
*              If the distance is less than aepsge the last vertex is 
*              saved in array.
*
* ASSUMTIONS  :kk multiple knots in start and end of parameterintervall.
*
* USE        :       
*
* REFERENCES :
*
*-                                                 
* CALLS      : s1840,s1710,s1612(recursive),s1235,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 31. Jan 1989
*
*********************************************************************
*/
{
  int kstat;          /* Status variable                                 */
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kdim;           /* The dimension of the space in which the curve
			 lies. Equivalently, the number of components
			 of each B-spline coefficient.                   */
  int ki;             /* Local counter                                   */
  int kleft;          /* Pointer into knot vector array                  */
  int knbreak=0;      /* No. of kk-1 multiple knots                      */
  int knbpnt;         /* Local for jnbpnt                                */
  int kleng;          /* Local for jleng                                 */
  int kvlast;         /* Position of last vertex in vertex array         */
  int kpos=0;         /* Position of error                               */
  
  double *spoint=SISL_NULL;/* Pointer to array of points                      */
  double tdist;       /* Distance                                        */
  double tpar;        /* A parameter value of the curve                  */
  double *sbreak = SISL_NULL;  /* Array containing kk-1 multiple knot         */
  double *st;         /* Pointer to the first element of the knot vector
			 of the curve. The knot vector has [kn+kk]
			 elements.                                       */
  
  
  SISLCurve *qcnew1=SISL_NULL; /* Pointer to first new  curve-object       */
  SISLCurve *qcnew2=SISL_NULL; /* Pointer to second new curve-objec        */
  
  
  /* Check input   */
  
  if (aepsge <= (double)0.0) goto err120;
  
  /* Make locals  */
  
  spoint = *gpoint;
  knbpnt = *jnbpnt;
  kleng  = *jleng;
  
  /* Describe curve with local parameters.  */
  
  kn    = pc -> in;
  kk    = pc -> ik;
  kdim  = pc -> idim;          
  st    = pc -> et;
  
  /* Find all kk-1 multiple knots including start and endpoint */
  
  s1235(st,kn,kk,&knbreak,&sbreak,&kstat);
  if (kstat < 0) goto error;
  
  /* Always split curve at kk-1 multiple knots */
  
  if (knbreak > 2) 
    {
      for (ki=1; ki<knbreak-1; ki++)
	{
	  tpar = sbreak[ki];
	  
	  s1710(pc,tpar,&qcnew1,&qcnew2,&kstat);
	  if(kstat < 0) goto error;
	  
	  /* Recursion on first part */
	  
	  if (qcnew1)
	  {
	     s1612(qcnew1,aepsge,&spoint,&knbpnt,&kleng,&kstat);
	     if (kstat < 0) goto error;
	  }
	  
	  /* Recursion on second part */
	  
	  if (qcnew2)
	  {
	     s1612(qcnew2,aepsge,&spoint,&knbpnt,&kleng,&kstat);
	     if (kstat < 0) goto error;
	  }
	}
    }
  else      
    /* If no internal kk-1 multiple knots */
    {
      /* Find distance between control polygon and a straight line beetween
	 first and last vertex*/
      
      s1840(pc,&tdist,&kstat);
      if (kstat < 0) goto error;
      
      if (tdist < aepsge)
	{
	  /* Save last vertex */
	  
	  kvlast = (kn-1) * kdim;
	  
	  knbpnt += 1;
	  
	  /* Allocate place for 100 more points in gpoint if not enought place */
	  
	  if (kleng < (knbpnt+1)*kdim )
	    {
	      kleng  += 100*kdim;
	      spoint = increasearray(spoint,kleng,DOUBLE);
	      if (!spoint) goto err101;
	    }
	  
	  /* Save point in array */
	  
	  memcopy (&spoint[(knbpnt-1)*kdim],&pc->ecoef[kvlast],kdim,DOUBLE);
	}
      else
	{
	  /* Split curve in two parts. If any internal knots, split in the 
	     midlemost */
	  
	  tpar = (st[0] + st[kn+kk-1]) / (double)2.0; 
	  
	  if (kn > kk) 
	    {    
	      /* Localize the parameter value */
	      
	      kleft = 0;
	      s1219(st,kk,kn,&kleft,tpar,&kstat);
	      if (kstat < 0) goto error;
	      
	      /* Find the knot nearest to the midle of the parameter interval */
	      
	      if (fabs(tpar-st[kleft]) < fabs(st[kleft+1]-tpar)) 
		{
		  tpar = st[kleft];
		}
	      else
		{
		  tpar = st[kleft+1];
		}
	    }
	  
	  s1710(pc,tpar,&qcnew1,&qcnew2,&kstat);
	  if(kstat < 0) goto error;
	  
	  /* Recursion on first part */
	  
	  if (qcnew1)
	  {
	     s1612(qcnew1,aepsge,&spoint,&knbpnt,&kleng,&kstat);
	     if (kstat < 0) goto error;
	  }
	  
	  /* Recursion on second part */
	  
	  if (qcnew2)
	  {
	     s1612(qcnew2,aepsge,&spoint,&knbpnt,&kleng,&kstat);
	     if (kstat < 0) goto error;
	  }
	}
    }
  
  *gpoint = spoint;
  *jnbpnt = knbpnt;
  *jleng  = kleng;
  *jstat = 0;
  goto out;
  
  /* Error in memory allocation */
  
 err101: 
  *jstat = -101;
  s6err("s1612",*jstat,kpos);
  goto out;
  
  /* Error in input. Relative tollerance <=0  */
  
 err120: 
  *jstat = -120;
  s6err("s1612",*jstat,kpos);
  goto out;
  
  
  /* Error in lower level function */  
  
 error:  
  *jstat = kstat;
  s6err("s1612",*jstat,kpos); 
  goto out;
  
  
 out:
  
  /* Free space occupied by local arrays and objects.  */
  
  if (sbreak) freearray(sbreak);
  if (qcnew1) freeCurve(qcnew1);
  if (qcnew2) freeCurve(qcnew2);
  
  return;
}       
   
          
