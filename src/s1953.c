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
 * $Id: s1953.c,v 1.3 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1953

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1953(SISLCurve *pcurve,double epoint[],int idim,double aepsco,double aepsge,
	   int *jpt,double **gpar,int *jcrv,SISLIntcurve ***wcurve,int *jstat)
#else
void s1953(pcurve,epoint,idim,aepsco,aepsge,jpt,gpar,jcrv,wcurve,jstat)
     SISLCurve    *pcurve;
     double   epoint[];
     int      idim;
     double   aepsco;
     double   aepsge;
     int      *jpt;
     double   **gpar;
     int      *jcrv;
     SISLIntcurve ***wcurve;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find the closest point between the curve pcurve and the
*              point epoint.
*
*
*
* INPUT      : pcurve - Pointer to the curve in the closest point problem.
*              epoint - The point in the closest point problem.
*              idim   - Dimension of the space in which epoint lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : jpt    - Number of single closest points.
*              gpar   - Array containing the parameter values of the
*                       single closest points in the parameter interval
*                       of the curve. The points lie continuous.
*                       Closest curves are stored in wcurve.
*              jcrv   - Number of closest curves.
*              wcurve - Array containing descriptions of the closest
*                       curves. The curves are only described by points
*                       in the parameter interval. The curve-pointers points
*                       to nothing. (See descrjption of Intcurve
*                       in intcurve.dcl).
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : Put the equation of the curve into the equation of a
*              circel/sphere with center epoint and radius zero. The
*              result is a curve of dimension one. Find the minimum
*              points of this curve.
*
*
*
* REFERENCES :
*
*-
* CALLS      : s1321,s1370,s1920,freeCurve.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, 1994-08.  Updated to handle
*              input curve with periodic basis correctly.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                    */
  int kpos = 0;             /* Position of error.                        */
  int kdim = 1;             /* Dimension of curve in extremal problem.   */
  double tradius = 0;       /* Radius of circel/sphere describing point. */
  double tdir = -1;         /* Direction of extremal value.              */
  double sarray[16];        /* Matrix describing circel/sphere.          */
  SISLCurve *qc = SISL_NULL;     /* Curve of which to find extremal points.   */
  SISLCurve *qkreg = SISL_NULL;  /* Input curve with ensured k-regular basis. */
  int ratflag = 0;          /* Flag to indicate if curve is rational.    */
  int ki;                   /* Counter.                                  */


  /* Ensure input has k-regular knot basis. */

  if ( pcurve -> cuopen == SISL_CRV_PERIODIC )
  {
    /* Cyclic (periodic) curve */

    make_cv_kreg(pcurve, &qkreg, &kstat);
    if ( kstat < 0 )  goto error;
  }
  else
    qkreg = pcurve;


  /* Test input.  */

  if ( idim != 2 && idim != 3 )  goto err105;
  if ( qkreg -> idim != idim )  goto err106;

  if ( qkreg -> ikind == 2 || qkreg -> ikind == 4 )  ratflag = 1;

  /* Make a matrix of dimension (idim+1)*(idim+1) to describe
     the circel/shpere.                                        */

  s1321(epoint, tradius, idim, kdim, sarray, &kstat);
  if ( kstat < 0 )  goto error;

  /* Put curve into circel/sphere equation.  */

  s1370(qkreg, sarray, idim, kdim, ratflag, &qc, &kstat);
  if ( kstat < 0 )  goto error;

  /* Find minimum points of the new curve.  */

  s1920(qc, &tdir, kdim, aepsco, aepsge, jpt, gpar, jcrv, wcurve, &kstat);
  if ( kstat < 0 )  goto error;

  /* Handle periodicity (remove extraneous points) */
  if ( *jpt > 1  &&  idim > 1  &&  pcurve -> cuopen == SISL_CRV_PERIODIC )
  {
    for ( ki=0;  ki < (*jpt);  ki++ )
    {
      if ( (*gpar)[ki] == pcurve -> et[pcurve->in] )
      {
	(*jpt)--;
	(*gpar)[ki] = (*gpar)[*jpt];
	ki--;
      }
    }
  }

  /* Closest points/intervals found.  */

  *jstat = 0;
  goto out;


  /* Error in input. Dimension not equal to 2 or 3.  */

 err105:
  *jstat = -105;
  s6err("s1953",*jstat,kpos);
  goto out;

  /* Dimensions conflicting.  */

 err106:
  *jstat = -106;
  s6err("s1953",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

 error:
  *jstat = kstat;
  s6err("s1953",*jstat,kpos);
  goto out;

 out:

  /* Free allocated space.  */

  if ( qkreg  &&  qkreg != pcurve )  freeCurve(qkreg);
  if (qc)  freeCurve(qc);

  return;
}
