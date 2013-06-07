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
 * $Id: s1513.c,v 1.2 2001-03-19 15:58:50 afr Exp $
 *
 */


#define S1513

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1513 (SISLSurf * psurf, double qpoint[], double bvec[], int idim, 
	SISLSurf ** rsurf, int *jstat)
#else
void 
s1513 (psurf, qpoint, bvec, idim, rsurf, jstat)
     SISLSurf *psurf;
     double qpoint[];
     double bvec[];
     int idim;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To make the function whose zeroes describes the circular silhouette
*              lines of a surface. The circular silhouette lines are defined
*              by a point qpoint and an axis bvec.
*                f(u,v) = N(u,v) x (P(u,v) - Q) . B = 0
*
*
*
* INPUT      : psurf  - Pointer to input surface
*              qpoint  - Point on axis vector.
*              bvec  - Spin axis vector.
B
*
*
*
* OUTPUT     : rsurf  - The resulting surface
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We first make the appropriate knot vector, then we calulate
*              parameter values for the interpolation, then the positions and normal
*              vectors of the surface are calculated at these parameter
*              values. The function f(U,V)
*              is calculated at the parameter values and these f values are interpolated
*              to put f, which describes the silhouette lines, into B-spline form.
*
*              For nonrational surfaces the function whose zeroes form the
*              silhouette is
*
*                      f = N(U,V) x (P(U,V) - Q) . B  = (P xP ) x (P - Q) . B
*                                                         U  V
*
*              which, if P is degree m x n, is degree 3m-1 x 3n-1
*              or order 3(m+1)-3 x 3(n+1)-3.
*              In the rational case, setting  P(U,V) = T(U,V) / w(U,V), we find
*
*                                                     4
*                      P xP =(wT -w T) x (wT -w T) / w
*                       U  V    U  U        V  V
*
*                                                      3
*                           =(wT xT +w T xT+w TxT ) / w
*                               U  V  U V    V   U
*
*              so we find the zeroes of
*
*                      f = (wT xT +w T xT+w TxT ) x (T - wQ) . B
*                             U  V  U V    V   U
*
*              which, if P is degree m x n, is degree 4m-1 x 4n-1
*              or order 4(m+1)-4 x 4(n+1)-4.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Mike Floater, SI, 1990-11
* REVISED BY : Mike Floater, SI, 1991-04 for rational surfaces
*
*********************************************************************
*/
{
  double cross[3];		/* Cross product of norm and diff                   */
  double diff[3];		/* Difference between surface point and eye point   */
  int kn1;			/* Number of vertices of psurf in first par.dir     */
  int kk1;			/* Order in  psurf in first par.dir                 */
  int kn2;			/* Number of vertices of psurf in second par.dir    */
  int kk2;			/* Order in  psurf in second par.dir                */
  int kjkk1;			/* Order of interpolated basis in first par.dir     */
  int kjkn1;			/* Number of vertices in interpolated basis first.dr*/
  int kjkk2;			/* Order of interpolated basis in first par dir     */
  int kjkn2;			/* Number of vertices in interpolated basis secnd.dr*/
  int kdim;			/* Number of dimesions in psurf                     */
  int kstat;			/* Local status variable                            */
  int kpos = 0;			/* Position indicator for errors                    */
  int kzero = 0;		/* Value 0                                          */
  int kone = 1;			/* Value 1                                          */
  int cuopen;			/* Open/Closed flag                                 */
  int ki, kj;	     	        /* Loop control variable                            */
  int kp;			/* Index of points put into conic equation          */
  int klfs = 0;			/* Pointer into knot vector                         */
  int klft = 0;			/* Pointer into knot vector                         */

  double *st1 = SISL_NULL;		/* First knot vector is psurf                       */
  double *st2 = SISL_NULL;		/* Second knot vector is psurf                      */
  double *sval1 = SISL_NULL;		/* Array of values of surface put into torus eq.    */
  double *sval2 = SISL_NULL;
  double *sval3 = SISL_NULL;

  double *sgt1 = SISL_NULL;		/* Knot vector in first parameter direction of
				   surface put into torus equation                  */
  double *sgt2 = SISL_NULL;		/* Knot vector in second parameter direction of
				   surface put into torus equation                  */
  double sder[12];		/* Point on the surface                             */
  double spar[2];		/* Current parameter pair                           */
  double snorm[3];		/* Normal vector                                    */
  int i;			/* a loop variable                                  */
  int ikind;			/* kind of surface psurf is                         */
  double tutv[3];		/* T_u x T_v  in rational case                      */
  double tvt[3];		/* T_v x T    in rational case                      */
  double ttu[3];		/* T   x T_u  in rational case                      */
  SISLSurf *tempsurf = SISL_NULL;	/* only used for rational surfaces                  */
  double *par1=SISL_NULL;
  double *par2=SISL_NULL;
  int    *der1=SISL_NULL;
  int    *der2=SISL_NULL;


  if (idim != psurf->idim) goto err104;

  /* Make local pointers */

  kn1 = psurf->in1;
  kk1 = psurf->ik1;
  kn2 = psurf->in2;
  kk2 = psurf->ik2;
  kdim = psurf->idim;
  st1 = psurf->et1;
  st2 = psurf->et2;
  ikind = psurf->ikind;

  /* Test input */

  if (kdim != 3)
    goto err104;

  if (ikind == 2 || ikind == 4)
    {
      /* A tricky way to evaluate the derivatives of the HOMOGENEOUS form
	 of psurf. In other words we need the derivs of T(u,v) where
	 p(u,v) = T(u,v) / w(u,v).
	 We should really have a separate evaluator for this
	 but I didn't want to mess around with the existing evaluator
	 which does the division automatically. MF 16/4/91 
     */
     
      tempsurf = newSurf (kn1, kn2, kk1, kk2, st1, st2, psurf->rcoef, ikind - 1, kdim + 1, 0);
      if (tempsurf == SISL_NULL) goto err171;
			    
      tempsurf->cuopen_1 = psurf->cuopen_1;			    
      tempsurf->cuopen_2 = psurf->cuopen_2;			    
      kjkk1 = 4 * kk1 - 4;
      kjkk2 = 4 * kk2 - 4;
    }
  else
    {
      tempsurf = psurf;
      kjkk1 = 3 * kk1 - 3;
      kjkk2 = 3 * kk2 - 3;
    }

  /* Make description of knot array for interpolation in first parameter
     direction */

  s1381 (st1, kn1, kk1, &sgt1, &kjkn1, kjkk1, &kstat);
  if (kstat < 0)
    goto error;

  /* Make parameter values and derivative indicators */

  s1890(sgt1,kjkk1,kjkn1,&par1,&der1,&kstat);
  if (kstat < 0)
    goto error;


  /* Make description of knot array for interpolation in second parameter
     direction
     */

  s1381 (st2, kn2, kk2, &sgt2, &kjkn2, kjkk2, &kstat);
  if (kstat < 0)
    goto error;

  /* Make parameter values and derivative indicators */

  s1890(sgt2,kjkk2,kjkn2,&par2,&der2,&kstat);
  if (kstat < 0)
    goto error;

  /* Allocate array for values of surface put into torus equation */

  sval1 = newarray (kjkn1 * kjkn2, DOUBLE);
  if (sval1 == SISL_NULL)
    goto err101;

  /* Calculate values to be interpolated */

  /* Index of point to be stored */

  kp = 0;

  for (kj = 0; kj < kjkn2; kj++)
    {

      spar[1] = par2[kj];

      for (ki = 0; ki < kjkn1; ki++)
	{
	  /* Calculate values on 3-D surface */

	  spar[0] = par1[ki];

	  s1421 (tempsurf, 1, spar, &klfs, &klft, sder, snorm, &kstat);
	  if (kstat < 0)
	    goto error;

	  if (ikind == 2 || ikind == 4)
	    {
	      /* Calculate (wT xT +w T xT+w TxT )
	                      U  V  U V    V   U
                 instead of normal to surface in rational case. */

	      s6crss (sder + 4, sder + 8, tutv);
	      s6crss (sder + 8, sder, tvt);
	      s6crss (sder, sder + 4, ttu);
	      for (i = 0; i < 3; i++)
		{
		  snorm[i] = sder[3] * tutv[i] + sder[7] * tvt[i] + sder[11] * ttu[i];
		  diff[i] = sder[i] - sder[3] * qpoint[i];
		}
	    }
	  else
	    {
	      s6diff (sder, qpoint, kdim, diff);
	    }

	  /* Make function f(u,v) = N(u,v) x ( P(u,v) - Q ) . B   */

	  s6crss (snorm, diff, cross);
	  sval1[kp++] = s6scpr (cross, bvec, kdim);
	}
    }

  cuopen = TRUE;

  /* Interpolate in second parameter direction, the first parameter direction
     is treated as a point of dimension kjkn1 */

  s1891(par2,sval1,kjkn1,kjkn2,kone,der2,cuopen,sgt2,&sval2,&kjkn2,kjkk2,kzero,kzero,&kstat);
  if (kstat < 0)
    goto error;


  /* Interpolate in first parameter direction, perform kjkn2 interpolations
     of one dimensional data */

  s1891(par1,sval2,kone,kjkn1,kjkn2,der1,cuopen,sgt1,&sval3,&kjkn1,kjkk1,kzero,kzero,&kstat);
  if (kstat < 0)
    goto error;

  if (kstat < 0)
    goto error;

  *rsurf = SISL_NULL;
  *rsurf = newSurf (kjkn1, kjkn2, kjkk1, kjkk2, sgt1, sgt2, sval3, 1, 1, 1);

  if (*rsurf == SISL_NULL)
    goto err171;

  (*rsurf)->cuopen_1 = psurf->cuopen_1;  
  (*rsurf)->cuopen_2 = psurf->cuopen_2; 
  
  /* Ok ! */


  *jstat = 0;
  goto out;

  /* Error in lower level function */

error:
  *jstat = kstat;
  s6err ("s1513", *jstat, kpos);
  goto out;

  /* Error in space allocation */
err101:
  *jstat = -101;
  s6err ("s1513", *jstat, kpos);
  goto out;

  /* Error in dimension */
err104:
  *jstat = -104;
  s6err ("s1513", *jstat, kpos);
  goto out;

  /* Could not create surfe */

err171:
  *jstat = -171;
  s6err ("s1513", *jstat, kpos);
  goto out;

out:

  /* Release allocated arrays */

  if (sgt1 != SISL_NULL)
    freearray (sgt1);
  if (sgt2 != SISL_NULL)
    freearray (sgt2);
  if (sval1 != SISL_NULL)
    freearray (sval1);
  if (sval2 != SISL_NULL)
    freearray (sval2);
  if (sval3 != SISL_NULL)
    freearray (sval3);
  if (par1 != SISL_NULL)
    freearray(par1);
  if (par2 != SISL_NULL)
    freearray(par2);
  if (der1 != SISL_NULL)
    freearray(der1);
  if (der2 != SISL_NULL)
    freearray(der2);
  if ((ikind == 2 || ikind == 4) && (tempsurf != SISL_NULL))
    freeSurf (tempsurf);

  return;
}
