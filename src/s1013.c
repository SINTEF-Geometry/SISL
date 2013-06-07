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
 * $Id: s1013.c,v 1.2 2001-03-19 15:58:40 afr Exp $
 *
 */


#define S1013

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s1013(SISLCurve *pcurve, double ang, double ang_tol,
	    double guess_par, double *iter_par, int *jstat)
#else
      void s1013(pcurve, ang, ang_tol, guess_par, iter_par, jstat)
	 SISLCurve   *pcurve;
	 double       ang;
	 double       ang_tol;
	 double       guess_par; 
	 double      *iter_par;
	 int         *jstat;
#endif
/********************************************************************
*                                                                   
* PURPOSE    : Find a point on a 2-Dimensional B-spline curve that
*              has a given direction.
*
*
*
* INPUT      : pcurve    - Pointer to the curve.
*              ang       - The angle (in radians) describing the wanted direction.
*              ang_tol   - The anular tolerance (in radians).
*              guess_par - Start parameter value on the B-spl crv.
*
*
* OUTPUT     : iter_par  - The found parameter value on the B-spl crv.
*              jstat     - status messages  
*                                = 2   : A minimum distanse found.
*                                = 1   : Intersection found.
*                                < 0   : error.
*
*
* METHOD     : Convert the problem to a onedimention zero problem
*              by making the first derivative of the curve and take
*              set the ratio of the components equal to the angle.
*              
*
*
* REFERENCES :
*
*- 
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SI, 92-03.
*
*********************************************************************
*/                                                               
{                                                                     
  int kstat   = 0;             /* Local status variable.                       */
  int kpos    = 0;             /* Position of error.                           */
  int kdim    = 2;             /* Legal dimension.                             */
  int ki      = 0;             /* Loop control.                                */
  int derive  = 1;             /* The first derivative.                        */  
  int kleft   = 0;             /* Knot navigator.                              */  
  double zero = 0;             /* Value to iteration.                          */ 
  double si,co;                /* Sin and cosin to the angle.                  */
  double *dim_one, *dim_two;   /* Pointers into the vertice array.             */
  double iter_ang;             /* The error angle                              */
  double sder[4];              /* Result from s1221.                           */
  double help_arr[2];          /* Help array containing angle as direction     */
  SISLCurve *testcurve = SISL_NULL; /* The curve to iterate on.                     */
  SISLPoint *p1 = newPoint (&zero, 1, 0);
  /* ------------------------------------------------------------------------- */
  
  /* Test input. ---------------------------------------------- */
  if (!p1) goto err101;
  if (pcurve->idim != kdim) goto err105;
  
  /* Find the ratio between the y and x direction. -------------*/
  help_arr[0] = co = cos(ang);
  help_arr[1] = si = sin(ang);
  
  /* Convert the problem. ------------------------------------- */
  s1720(pcurve, derive, &testcurve, &kstat);
  if (kstat < 0) goto error;
  
  testcurve->idim--;
  dim_one = testcurve->ecoef;
  dim_two = testcurve->ecoef;
  for (ki = 0; ki < testcurve->in; ki++)
    {
       *dim_one = co*(*(dim_two + 1)) - si*(*dim_two);
       dim_one++;
       dim_two += 2;
    }
  
  
  /* Iteration. ---------------------------------------------- */
  s1771(p1, testcurve, REL_COMP_RES, testcurve->et[testcurve->ik - 1],
	testcurve->et[testcurve->in], guess_par, iter_par, &kstat);
  if (kstat < 0) goto error;
  
  /* Analyse the result */
  s1221(pcurve, derive=1, *iter_par, &kleft, sder, &kstat);
  if (kstat < 0) goto error;

  iter_ang = s6ang(sder+2, help_arr, kdim);
  
  if (iter_ang < ang_tol)   *jstat = 1;
  else   *jstat = 2;

  goto out;
  

  /* EXITS. ------------------------------------------------ */

  /* Error in space allocation.  */
 err101: *jstat = -101;
  s6err("s1013",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2.  */
 err105: *jstat = -105;
  s6err("s1013",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  error : *jstat = kstat;
  s6err("s1013",*jstat,kpos);
  goto out;
  
 out:
  /* Free allocated space.  */
  if (testcurve) freeCurve(testcurve);
  if (p1)        freePoint(p1);
 
  return;
}                                               
                                           
                       
