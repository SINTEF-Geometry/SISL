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
 * $Id: s1226.c,v 1.4 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1226

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1226(SISLCurve *curve,int der,double parvalue,int *leftknot,
      double derive[],double curvature[], double *radius_of_curvature,
      int *jstat)
#else
void s1226(curve,der,parvalue,leftknot,derive,curvature,
	   radius_of_curvature,jstat)
     SISLCurve  *curve;
     int    der;
     double parvalue;
     int    *leftknot;
     double derive[];
     double curvature[];
     double *radius_of_curvature;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Evaluate position, first derivative, curvature and radius of
*              curvature of a curve at a given parameter value, from the
*              right hand side.
*
*
*
* INPUT      : curve    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              der      - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              parvalue - The parameter value at which to compute
*                       position and derivatives.
*
*
*
* INPUT/OUTPUT : leftknot - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*
*                          et[ileft] < parvalue <= et[ileft+1]
*
*                        should hold. (If parvalue == et[ik-1] then ileft
*                        should be ik-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : derive   - Double array of dimension [(ider+1)*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*                       (The C declaration of eder as a two dimensional array
*                       would therefore be eder[ider+1,idim].)
*              curvature - Array of dimension idim
*              radius_of_curvature -  A radius of curvature =-1, indicates
*                       that the radius of curvature is infinit.
*              jstat  - Status messages
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : The derivatives are evaluated from the right hand
*              side by  s1221
*              The curvature and the radius of curvature are evaluated
*              by s1307
*
* REFERENCES :
*
*-
* CALLS      : s1221, s1307
*
* WRITTEN BY : Cathrine Tegnander
* MODIFIED BY :
* REVISED BY :
*
*********************************************************************
*/
{
  int ider=2;                 /* minimum number of derivatives needed to
                                 calculated the curvature */
  int kdim = curve -> idim;   /* copy curve attribute to local parameter */
  int kstat=0;                /* local status variable */
  int kpos = 0;               /* local error position  */
  int iknot = 0;              /* local version of leftknot */
  double *iderive = SISL_NULL;     /* pointer to array used to store the position
                                 and the first and second derivatives  */
  double *egeo = SISL_NULL;        /* pointer to store curvature and radius of
                                 curvature */

  iderive= newarray(3*kdim,DOUBLE);
  if (iderive == SISL_NULL)
    goto err101;
  egeo= newarray(3*kdim+1,DOUBLE);
  if (egeo == SISL_NULL)
    goto err101;


  /* Evaluate the derivatives */
    if (der<2)
      {
	s1221(curve,ider,parvalue,&iknot,iderive,&kstat);
	if (kstat<0) goto error;

	/* Copy position */
	memcopy(derive,iderive,(der+1)*kdim,DOUBLE);
      }
  else
    {
      s1221(curve,der,parvalue,&iknot,derive,&kstat);
      if (kstat<0) goto error;

      	/* Copy position */
	memcopy(iderive,derive,3*kdim,DOUBLE);
    }
    *leftknot = iknot;

  /* Evaluate the curvature and the radius_of_curvature */
  s1307(iderive,kdim,egeo,&kstat);
  if (kstat<0) goto error;

  /* Copy position */
  memcopy(curvature,egeo+kdim*2,kdim,DOUBLE);

  *radius_of_curvature = egeo[kdim*3];

  /* Free memory */
 freearray(iderive);
 freearray(egeo);

 /* Successful computations.  */

  *jstat = 0;
  goto out;

  /* Error in scratch allocation */
 err101:
  *jstat = -101;
  s6err("s1226",*jstat,kpos);

  /* Error in lower level routine.  */

 error:  *jstat = kstat;
  s6err("S1226",*jstat,kpos);
  goto out;

 out: return;

}
