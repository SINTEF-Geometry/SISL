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
 * $Id: s6addcurve.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6ADDCURVE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static void s6addcurve_s9moveknots(double [],int,double,double,
				   double [],int *); 
#else
static void s6addcurve_s9moveknots(); 
#endif

#if defined(SISLNEEDPROTOTYPES)
void s6addcurve(SISLCurve *pc1,SISLCurve *pc2,int isign,
		SISLCurve **rcurve,int *jstat)
#else
void s6addcurve(pc1,pc2,isign,rcurve,jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     int isign;
     SISLCurve **rcurve;
     int *jstat;
#endif   
/*
*********************************************************************
* 
* PURPOSE    : Compute the sum or difference between two B-spline
*              curves depending on the constang isign. That is find
*              the curve pc1 + isign*pc2.
* 
* 
* 
* INPUT      : pc1      - First curve in expression.
*              pc2      - Second curve in expression.
*              isign    - Sign of second curve. 
*
* 
* OUTPUT     : rcurve   - Sum/difference between the input curves.
*              jstat    - status messages 
*                         > 0 : warning
*                         = 0 : ok 
*                         < 0 : error 
* 
* 
* METHOD     : Express the curves on the same knot vector. Compute the
*              sum difference between corresponding vertices.
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : s1932  - Express curve-set in given basis.  
*              s1363  - Pick parametrisation of curve.
*              s6takeunion - Union of two ordered vectors. 
*              make_cv_kreg - Make curve k-regular.
*              s1333_count - Count continuity in ends of closed curves.
*              make_cv_cyclic - Represent curve in a periodic basis.
*              newCurve  -   Create new curve-object. 
*              s6addcurve_s9moveknots() - Move knotvector into new parameter interval.
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 08.90.
* REWISED BY : Vibeke Skytt, SI, 05.92. To treat periodic curves.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Status variable.  */
  int ki;                /* Counter.          */
  int kdim = pc1->idim;  /* Dimension of geometry space.              */
  int knbcrv = 2;        /* Number of input curves.                   */
  int kknot1 = pc1->in + pc1->ik;  /* Number of knots of first curve. */
  int kknot2 = pc2->in + pc2->ik;  /* Number of knots of second curve. */
  int korder;            /* Order of output curve.     */
  int ktau;              /* Number of knots / number of vertices of 
			    output curve.                             */
  int kcont;             /* Continuity at ends in periodic case.      */
  int kopen;             /* Open/closed/periodic parameter.           */
  int kkind = pc1->ikind;  /* Kind of curves.                         */
  int kcopy = 1;         /* Copy arrays when creating new curve.      */
  double tmin1,tmax1;    /* Parameter interval of first input curve.  */
  double tmin2,tmax2;    /* Parameter interval of second input curve. */
  double *st = SISL_NULL;     /* Knot vector of second curve.              */
  double *stau = SISL_NULL;   /* Knot vector of output curve.              */
  double *scoef = SISL_NULL;  /* Vertices of input curves represented in 
			    the same knot vector, then of output curve. */
  SISLCurve *ucurves[2];     /* Local pointer array to input curves.      */
  
  ucurves[0] = ucurves[1] = SISL_NULL;
  
  /* Test input.  */

  if (kdim != pc2->idim) goto err106;
  
  /* Make curves k-regular.  */

  make_cv_kreg(pc1,ucurves,&kstat);
  if (kstat < 0) goto error;
 
  make_cv_kreg(pc2,ucurves+1,&kstat);
  if (kstat < 0) goto error;
  
  /* Allocate scratch for local knot vector.  */

  if ((st = newarray(kknot2,DOUBLE)) == SISL_NULL) goto err101;
  
  /* Fetch endparameter values of the input curves.  */
   
  s1363(ucurves[0],&tmin1,&tmax1,&kstat);
  if (kstat < 0) goto error;
  
  s1363(ucurves[1],&tmin2,&tmax2,&kstat);
  if (kstat < 0) goto error;
  
  if (DNEQUAL(tmin1,tmin2) || DNEQUAL(tmax1,tmax2))
    {
      /* Move the knot vector of the second curve into the parameter
	 interval of the first curve.  */

       s6addcurve_s9moveknots(ucurves[1]->et,kknot2,tmin1,tmax1,st,&kstat);
      if (kstat < 0) goto error;
    }
  else memcopy(st,ucurves[1]->et,kknot2,DOUBLE);
  
  /* Find the union of the knot vectors of the two curves.  */

  s6takeunion(ucurves[0]->et,kknot1,st,kknot2,&stau,&ktau,&kstat);
  if (kstat < 0) goto error;
  
  korder = MAX(ucurves[0]->ik,ucurves[1]->ik);
  ktau -= korder;
  
  /* Express the curves in the basis.  */

  s1932(2,ucurves,tmin1,tmax1,stau,ktau,korder,&scoef,&kstat);
  if (kstat < 0) goto error;
  
  /* Find the coefficients of the sum/difference curve.  */

  for (ki=0; ki<ktau*kdim; ki++)
    scoef[ki] += (double)isign*scoef[ktau*kdim+ki];
  
  /* Present the sum/difference as a curve.  */

  if ((*rcurve = newCurve(ktau,korder,stau,scoef,kkind,kdim,kcopy)) == SISL_NULL)
    goto err101;
  
  (*rcurve)->cuopen = kopen = MAX(pc1->cuopen,pc2->cuopen);
  if (kopen == SISL_CRV_PERIODIC)
  {
     /* Represent the output curve cyclic. First find continuity. */
     
     s1333_count(knbcrv,ucurves,&kcont,&kstat);
     if (kstat < 0) goto error;
     
     make_cv_cyclic(*rcurve,kcont,&kstat);
     if (kstat < 0) goto error;
  }
  
  /* Task performed.  */

  *jstat = 0;
  goto out;
  
  /* Error in scratch allocation.  */

  err101 :
    *jstat = -101;
  goto out;
  
  /* Error in input. Conflicting dimensions.  */

  err106 :
    *jstat = -106;
  goto out;
  
  /* Error in lower order routine.  */

  error :
    *jstat = kstat;
  goto out;
  
  out :

    /* Free scratch occupied by local arrays.  */

    if (st != SISL_NULL) freearray(st);
  if (stau != SISL_NULL) freearray(stau);
  if (scoef != SISL_NULL) freearray(scoef);
  if (ucurves[0] != SISL_NULL) freeCurve(ucurves[0]);		     
  if (ucurves[1] != SISL_NULL) freeCurve(ucurves[1]);		     
  
  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  s6addcurve_s9moveknots(double et[],int inmbknot,double astart,
			 double aend,double etmoved[],int *jstat)
#else
static void s6addcurve_s9moveknots(et,inmbknot,astart,aend,etmoved,jstat)
     double et[];
     int inmbknot;
     double astart;
     double aend;
     double etmoved[];
     int *jstat;
#endif     
/*
*********************************************************************
* 
* PURPOSE    : Move knot vector into the parameter interval
*              [astart,aend].  
* 
* 
* 
* INPUT      : et       - Knot vector.
*              inmbknot - Number of knots in the vector et.
*              astart   - Start parameter value of new parameter interval.
*              aend     - End parameter value of new parameter interval. 
*
* 
* OUTPUT     : etmoved  - Knot vector on new parameter interval.
*              jstat    - status messages 
*                         = 1 : Parameter interval changed. Thus,
*                               reparametrization is performed.
*                         = 0 : ok 
*                         < 0 : error 
* 
* 
* METHOD     : 
*
*
* REFERENCES : 
*              
*
* USE        :
*
*-
* CALLS      : 
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 07.90.
*
*********************************************************************
*/
{
  int kwarn = 0;             /* Indicates if reparametrization is performed. */
  double tstart = et[0];     /* Start parameter of original knot vector.     */
  double tend = et[inmbknot-1];   /* End parameter of original knot vector.  */
  double tfac = (aend - astart)/(tend - tstart);  /* Factor used when moving
						     knot vector.            */
  double *s1,*s2,*s3;        /* Pointers used to traverse knot vectors.      */
  
  /* Test if reparametrization is performed.  */

  if (DNEQUAL(aend-astart,tend-tstart)) kwarn = 1;
  
  /* Put knot vector into new parameter interval.  */

  for (s1=et,s2=et+inmbknot,s3=etmoved; s1<s2; s1++,s3++)
    *s3 = astart + (*s1 - tstart)*tfac;
  
  *jstat = kwarn;
  return;
}


