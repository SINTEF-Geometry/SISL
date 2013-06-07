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
 * $Id: s1604.c,v 1.3 2001-03-19 15:58:51 afr Exp $
 *
 */


#define S1604

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1604(double epoint[],int inbpnt,double astpar,int iopen,int idim,int ik,
	   SISLCurve **rc,int *jstat)
#else
void s1604(epoint,inbpnt,astpar,iopen,idim,ik,rc,jstat)
	   double epoint[];
	   int    inbpnt;
	   double astpar;
	   int    iopen;
	   int    idim;
	   int    ik;
	   SISLCurve  **rc;
	   int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To calculate a B-spline curve using the input points as
*              controlling vertices. The distances between the points are
*              used as parametrization.
*
*
* INPUT      : epoint - The array containing the points to be used as
*                       controlling vertices of the B-spline curve.
*              inbpnt - No. of points in epoint.
*              astpar - Parameter value to be used at the start of the curve.
*              iopen  - Open/close condition (Open=1,Close=0)
*              idim   - The dimension of the space
*              ik     - The order of the B-spline curve to be produced.
*
* OUTPUT     : rc     - Pointer to the curve
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : First the parametrization of the curve is calculated.
*              If more than ik adjacent vertices are equal then the
*              superfluous vertices are removed. Then the knots are
*              calculated.
*
* EXAMPLE OF USE:
*
* REFERENCES :
*
*-
* CALLS      : s6dist,s1902,s1713,s1750,newCurve,s6err
*
*
* WRITTEN BY : Qyvind Hjelle, SI, Oslo, Norway. 22. Nov 1988
* REVISED BY : Bjoern Olav Hoset, SI, Oslo, Norway, Feb. 1992
*              Calls s1902 instead of fortran functions.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway. 24/08-1994. Changed
*              icopy param. (from 0 to 2) to fix memory leak.
*********************************************************************
*/
{
  int kstat;          /* Status variable                                 */
  int kn;             /* The number of B-splines, i.e., the dimension of
			 the spline space associated with the knot
			 vector.                                         */
  int kk;             /* The polynomial order of the curve.              */
  int kpos=0;         /* Position of error                               */
  int ki;             /* Counter for loop control                        */

  double *spara=SISL_NULL; /* Pointer to parameterization array               */
  double *scoef=SISL_NULL; /* Pointer to vertex array                         */
  double *sknot=SISL_NULL; /* Pointer to knot vector                          */
  double tdist;       /* Distance */
  double tlastpar;    /* Last value in the parameterization array        */
  SISLCurve *qc=SISL_NULL;
  SISLCurve *qc2 = SISL_NULL;


  /* Control input     */

/* -> new statement guen & ujk Thu Jul  2 14:59:05 MESZ 1992 */

/* make compatible to old use of s1604 */
  if (iopen==SISL_CRV_CLOSED) iopen = SISL_CRV_PERIODIC;

/* a new version with input-iopen == rc->cuopen
   should be made, with a name different from any
   other.
  NOTE: There is an error in this function when iopen = 0
	qc as input to s1713 (and s1750) then has wrong flag !!*/

/* <- new statement guen & ujk Thu Jul  2 14:59:05 MESZ 1992 */

  kk = ik;
  if  (inbpnt < kk)
    kk = inbpnt;

  if  (kk < 2) goto err109;

  if  (iopen != SISL_CRV_OPEN && iopen != SISL_CRV_CLOSED
       && iopen != SISL_CRV_PERIODIC) goto err113;

  /* Allocate space for parameterization    */

  spara = newarray(inbpnt+1,DOUBLE);
  if (spara == SISL_NULL) goto err101;

  /* Calculate parameterization  */

  spara[0] = astpar;
  tlastpar = astpar;

  for (ki=1; ki<inbpnt; ki++)
    {
      tdist = s6dist(&epoint[ki*idim-idim],&epoint[ki*idim],idim);
      tlastpar = tlastpar + tdist;
      spara[ki] = tlastpar;
    }


  /* Calculate distance from first to last point and update
     parameterization array. To be useded if closed curve   */

  tdist = s6dist(epoint,&epoint[(inbpnt-1)*idim],idim);
  tlastpar = tlastpar + tdist;
  spara[inbpnt] = tlastpar;


  /* Find the knot vector     */

  s1902(spara,inbpnt,kk,iopen,&sknot,&kstat);
  if (kstat < 0 || sknot == SISL_NULL) goto error;

  /* Allocate space for verice array   */

  scoef = newarray((inbpnt+kk-1)*idim,DOUBLE);
  if (scoef == SISL_NULL) goto err101;

  /* Copy vertices */

  memcopy (scoef,epoint,inbpnt*idim,DOUBLE);
  kn = inbpnt;

  /* In case of closed curve, add the (kk-1) first points to the vertice
     array.    */

  if (!(iopen == SISL_CRV_OPEN))
    {
      memcopy(&scoef[inbpnt*idim],epoint,(kk-1)*idim,DOUBLE);
      kn = kn + kk - 1;
    }

  /* Make curve */
  /* VSK, MESZ. Do not copy arrays.
  qc = newCurve(kn,kk,sknot,scoef,1,idim,1);  */
  qc = newCurve(kn,kk,sknot,scoef,1,idim,2); /* icopy=2, PFU 24/08-94 */
  if (!qc) goto err101;

  qc->cuopen = iopen;

  if (iopen == SISL_CRV_CLOSED)
    {
      /* A closed, non-periodic curve is expected. Pick the part of the
	 interpolation curve that has got a full basis.  */

      s1713 (qc, sknot[kk - 1], sknot[kn], &qc2, &kstat);
      if (kstat < 0)
	goto error;

      if (qc != SISL_NULL)
	freeCurve (qc);
      qc = qc2;
    }

  /* Increase the order if the order was lowered when controlling input */

  if (kk < ik)
    {
				/* -> guen: inserted acc. to SCCS */
      s1750(qc,ik,rc,&kstat);
				/* <- guen: inserted acc. to SCCS */

				/* -> guen: removed acc. to SCCS  */
				/*      s1750(qc,ik,&qc,&kstat);  */
				/* <- guen: removed acc. to SCCS  */
      if (kstat< 0) goto error;
    }
				/* -> guen: inserted acc. to SCCS */
  else
    {
      *rc = qc;
      qc = SISL_NULL;
    }
				/* <- guen: inserted acc. to SCCS */
				/* -> guen: removed acc. to SCCS  */
				/*  if (qc) *rc = qc;             */
				/* -> guen: inserted acc. to SCCS */

  *jstat = 0;
  goto out;

  /* Error in memory allocation */

 err101:
  *jstat = -101;
  s6err("s1604",*jstat,kpos);
  goto out;

  /* Error in input, order less than 2 */

 err109:
  *jstat = -109;
  s6err("s1604",*jstat,kpos);
  goto out;

  /* Error in input, unknown kind of curve */

 err113:
  *jstat = -113;
  s6err("s1604",*jstat,kpos);
  goto out;


  /* Error in lower level function */

 error:
  *jstat = kstat;
  s6err("s1604",*jstat,kpos);
  goto out;

 out:
  if (spara  != SISL_NULL) freearray(spara);
  /* if (scoef  != SISL_NULL) freearray(scoef); (Freed by freeCurve(qc). */
  if (qc != SISL_NULL) freeCurve(qc);
  return;
}
