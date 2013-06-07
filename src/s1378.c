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
 * $Id: s1378.c,v 1.3 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1378

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1378 (SISLSurf * psurf, double econic[], int ideg, int idim,
       SISLSurf ** rsurf, int *jstat)
#else
void
s1378 (psurf, econic, ideg, idim, rsurf, jstat)
     SISLSurf *psurf;
     double econic[];
     int ideg;
     int idim;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To put a surface description into the descripiton of
*              a torus surface described by the input array econic.
*
* INPUT      : psurf  - Pointer to input surface
*              econic - Description of torus
*              ideg   - Type of conic: torus: ideg=1001
*              idim   - Dimension of object space
*
* OUTPUT     : rsurf  - The resulting surface
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : We first make the appropriate knot vector, then we calulate
*              parametervalues for the interpolation, then the appropriate
*              values of the surface put into the conic equation are found,
*              and at last the surface is interpolated.
*
* REFERENCES :
*
* CALLS      : s1376,s1890,s1424,s6scpr,s1891,s6err.
*
* WRITTEN BY : Tor Dokken, SI, 1988-11
* REVISED BY : Mike Floater, SI, 1991-01 for a rational surface.
*
*********************************************************************
*/
{
  int ikind;			/* type of surface psurf is                         */
  int kn1;			/* Number of vertices of psurf in first par.dir     */
  int kk1;			/* Order in  psurf in first par.dir                 */
  int kn2;			/* Number of vertices of psurf in second par.dir    */
  int kk2;			/* Order in  psurf in second par.dir                */
  int kjkk1;			/* Order of interpolated basis in first par.dir     */
  int kjkn1;			/* Number of vertices in interpolated basis first.dr*/
  int kjkk2;			/* Order of interpolated basis in first par SISLdir     */
  int kjkn2;			/* Number of vertices in interpolated basis secnd.dr*/
  int kdim;			/* Number of dimesions in psurf                     */
  int kstat;			/* Local status variable                            */
  int kpos = 0;			/* Position indicator for errors                    */
  int kzero = 0;		/* Value 0 needed in call s1891		          */
  int kone = 1;			/* Value 1 needed in call s1891			  */
  int cuopen;			/* Open/Closed flag                                 */
  int ki, kj, kl;		/* Loop control variable                            */
  int kp;			/* Index of points put into conic equation          */
  int klfs = 0;			/* Pointer into knot vector                         */
  int klft = 0;			/* Pointer into knot vector                         */
  double *st1 = SISL_NULL;		/* First knot vector is psurf                       */
  double *st2 = SISL_NULL;		/* Second knot vector is psurf                      */
  double *scentr = econic;	/* Center of torus             */
  double *saxis = econic + 3;	/* Axis of torus               */
  double tbigr = *(econic + 6);	/* Big radius of torus         */
  double tsmalr = *(econic + 7);/* Small radius of torus       */
  double tbigr2 = tbigr * tbigr;/* Square of big radius        */
  double tdiffr2 = tbigr2 - tsmalr * tsmalr;	/* Difference of square of radia*/
  double *sval1 = SISL_NULL;		/* Array of values of surface put into torus eq.    */
  double *sval2 = SISL_NULL;
  double *sval3 = SISL_NULL;
  double *sgt1 = SISL_NULL;		/* Knot vector in first parameter direction of
				   surface put into torus equation                  */
  double *sgt2 = SISL_NULL;		/* Knot vector in second parameter direction of
				   surface put into torus equation                  */
  double sy[3];			/* Difference between point and torus center        */
  double tzn;			/* Projection of sy onto torus axis                 */
  double tyy;			/* Square of length of sy                           */
  double tzz;			/* Square of length of sz                           */
  double ty;			/* Component of sy                                  */
  double tz;			/* Component of sz                                  */
  double sder[4];		/* SISLPoint on the surface                         */
  double spar[2];		/* Current parameter pair                           */
  double ww;			/* the weight of sder squared if psurf is rational  */
  double *par1 = SISL_NULL;		/* Parameter vaues in direction 1. 		  */
  double *par2 = SISL_NULL;		/* Parameter vaues in direction 2. 		  */
  int *der1 = SISL_NULL;		/* Derivative indicators in direction 1.		  */
  int *der2 = SISL_NULL;		/* Derivative indicators in direction 2.		  */
  SISLSurf *tempsurf = SISL_NULL;	/* only used for rational surfaces             */

  *jstat = 0;


  /* Test if torus. */

  if (ideg != 1001)
    goto err180;

  if (idim != psurf->idim)
    goto err104;

  /* Make local pointers. */

  kn1 = psurf->in1;
  kk1 = psurf->ik1;
  kn2 = psurf->in2;
  kk2 = psurf->ik2;
  kdim = psurf->idim;
  st1 = psurf->et1;
  st2 = psurf->et2;
  ikind = psurf->ikind;

  if (ikind == 2 || ikind == 4)
    {
      tempsurf = newSurf (kn1, kn2, kk1, kk2, st1, st2,
			  psurf->rcoef, ikind - 1, kdim + 1, 0);
      if (tempsurf == SISL_NULL)
	goto err171;
      tempsurf->cuopen_1 = psurf->cuopen_1;
      tempsurf->cuopen_2 = psurf->cuopen_2;
    }
  else
    {
      tempsurf = psurf;
    }

  /* Test input. */

  if (kdim != 3)
    goto err104;


  /* Make description of knot array for interpolation in first parameter
     direction. */

  s1376 (st1, kn1, kk1, &sgt1, &kjkn1, &kjkk1, &kstat);
  if (kstat < 0)
    goto error;


  /* Make parameter values and derivative indicators. */

  s1890 (sgt1, kjkk1, kjkn1, &par1, &der1, &kstat);
  if (kstat < 0)
    goto error;


  /* Make description of knot array for interpolation in second parameter
     direction. */

  s1376 (st2, kn2, kk2, &sgt2, &kjkn2, &kjkk2, &kstat);
  if (kstat < 0)
    goto error;


  /* Make parameter values and derivative indicators. */

  s1890 (sgt2, kjkk2, kjkn2, &par2, &der2, &kstat);
  if (kstat < 0)
    goto error;


  /* Allocate array for values of surface put into torus equation. */

  sval1 = newarray (kjkn1 * kjkn2, DOUBLE);
  if (sval1 == SISL_NULL)
    goto err101;


  /* Calculate values to be interpolated. */

  /* Index of point to be stored. */

  kp = 0;

  for (kj = 0; kj < kjkn2; kj++)
    {

      spar[1] = par2[kj];

      for (ki = 0; ki < kjkn1; ki++)
	{
	  /*  Calculate values on 3-D surface */

	  spar[0] = par1[ki];

	  s1424 (tempsurf, 0, 0, spar, &klfs, &klft, sder, &kstat);
	  if (kstat < 0)
	    goto error;

	  /*
	   *       The calculation of a point on the torus surface
	   *		 can be done in the following way.
	   *
	   *          y = p - scentr
	   *          z = y - (y saxis) saxis
	   *
	   *       The equation of the torus can be written
	   *
	   *                              2    2
	   *          (y - R z/sqrt(z z) )  - r = 0
	   *
	   *
	   *       or by elliminating the square root:
           *
           *          f =
	   *
	   *              2           2  2      2       2  2 2
	   *          (yy)  + 2 (yy)(R -r ) - 4R zz + (R -r )  = 0
	   *
           *       or in 4-D homogeneous coordinates:
           *
           *                                               4
           *          f =
	   *
	   *              2      2      2  2      2 2      4  2  2 2
	   *          (yy)  + 2 w (yy)(R -r ) - 4w R zz + w (R -r )  = 0
	   *
	   *         where Y = T - w*scentr,  p+T/w
	   *
	   *       We thus need to calculate yy and zz:
	   */

	  if (ikind == 2 || ikind == 4)
	    {
	      for (kl = 0; kl < 3; kl++)
		sy[kl] = sder[kl] - sder[3] * scentr[kl];
	      ww = sder[3] * sder[3];
	    }
	  else
	    {
	      for (kl = 0; kl < 3; kl++)
		sy[kl] = sder[kl] - scentr[kl];
	      ww = (double) 1.0;
	    }

	  tzn = s6scpr (sy, saxis, 3);

	  tyy = (double) 0.0;
	  tzz = (double) 0.0;

	  /*      Make z and necessary derivatives of z */

	  for (kl = 0; kl < 3; kl++)
	    {
	      ty = sy[kl];
	      tz = ty - tzn * saxis[kl];
	      tyy += ty * ty;
	      tzz += tz * tz;
	    }

	  /*                                      2            2   2
	     Now tyy = yy and tzz = zz, tbigr2 = R ,tdiffr2 = R - r   */

	  sval1[kp++] = tyy * tyy + ((double) 2.0 * ww * tyy + ww * ww * tdiffr2) * tdiffr2
	    - (double) 4.0 *ww * tbigr2 * tzz;
	}
    }

  cuopen = TRUE;

  /* Interpolate in second parameter direction, the first parameter direction
     is treated as a point of dimension kjkn1 */

  s1891 (par2, sval1, kjkn1, kjkn2, kone, der2, cuopen, sgt2, &sval2,
	 &kjkn2, kjkk2, kzero, kzero, &kstat);
  if (kstat < 0)
    goto error;


  /* Interpolate in first parameter direction, perform kjkn2 interpolations
     of one dimensional data */

  s1891 (par1, sval2, kone, kjkn1, kjkn2, der1, cuopen, sgt1, &sval3,
	 &kjkn1, kjkk1, kzero, kzero, &kstat);
  if (kstat < 0)
    goto error;

  *rsurf = SISL_NULL;
  *rsurf = newSurf (kjkn1, kjkn2, kjkk1, kjkk2, sgt1, sgt2, sval3, 1, 1, 1);
  if (*rsurf == SISL_NULL)
    goto err171;
  (*rsurf)->cuopen_1 = psurf->cuopen_1;
  (*rsurf)->cuopen_2 = psurf->cuopen_2;

  /* Ok ! */

  goto out;


  /* Error in lower level function */

error:
  *jstat = kstat;
  s6err ("s1378", *jstat, kpos);
  goto out;

  /* Error in space allocation */

err101:
  *jstat = -101;
  s6err ("s1378", *jstat, kpos);
  goto out;

  /* Dimension not equal to 3  or confliciting dim  */

err104:
  *jstat = -104;
  s6err ("s1378", *jstat, kpos);
  goto out;

  /* Could not create surface. */

err171:
  *jstat = -171;
  s6err ("s1378", *jstat, kpos);
  goto out;

  /* Error in torus description */

err180:
  *jstat = -180;
  s6err ("s1378", *jstat, kpos);
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
