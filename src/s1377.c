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
 * $Id: s1377.c,v 1.3 2001-03-19 15:58:48 afr Exp $
 *
 */


#define S1377

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1377 (SISLCurve * pcurv, double econic[], int ideg, int idim,
       SISLCurve ** rcurv, int *jstat)
#else
void
s1377 (pcurv, econic, ideg, idim, rcurv, jstat)
     SISLCurve *pcurv;
     double econic[];
     int ideg;
     int idim;
     SISLCurve **rcurv;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To put a curve description into the descripiton of
*              a torus surface described by the input array econic.
*
*
* INPUT      : pcurv  - Pointer to input curve
*              econic - Description of torus
*              ideg   - Type of conic: torus: ideg=1001
*              idim   - Dimension of object space
*
*
* OUTPUT     : rcurv  - The resulting curve
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We first make the appropriate knot vector, then we calulate
*              parametervalues for the interpolation, then the appropriate
*              values of the curve put into the conic equation are found,
*              and at last the curve is interpolated.
*
* REFERENCES :
*
* CALLS      : s1376,s1890,s1221,s6scpr,s1891,s6err.
*
* WRITTEN BY : Tor Dokken, SI, 1988-11
* REVISED BY : Mike Floater, SI, 1991-04 for a rational curve.
*
*********************************************************************
*/
{
  int ikind;			/* Type of curve pcurv is                           */
  int kn;			/* Number of vertices of pcurv                      */
  int kk;			/* Order in  pcurv                                  */
  int kjkk;			/* Order of interpolated basis                      */
  int kjkn;			/* Number of vertices in interpolated basis         */
  int kdim;			/* Number of dimensions in pcurv                    */
  int kstat;			/* Local status variable                            */
  int kpos = 0;			/* Position indicator for errors                    */
  int kzero = 0;		/* Value 0 needed in call s1891    	          */
  int kone = 1;			/* Value 1 needed in call s1891		          */
  int cuopen;			/* Open/Closed flag                                 */
  int ki, kj;			/* Loop control variable                            */
  int kleft = 0;		/* Pointer into knot vector                         */
  int *der = SISL_NULL;		/* Derivate indicators. */
  double *st = SISL_NULL;		/* First knot vector is pcurv                       */
  double *scentr = econic;	/* Center of torus             */
  double *saxis = econic + 3;	/* Axis of torus               */
  double tbigr = *(econic + 6);	/* Big radius of torus         */
  double tsmalr = *(econic + 7);/* Small radius of torus       */
  double tbigr2 = tbigr * tbigr;/* Square of big radius        */
  double tdiffr2 = tbigr2 - tsmalr * tsmalr;	/* Difference of square of radii*/
  double *sval1 = SISL_NULL;		/* Array of values of curve put into torus eq.      */
  double *sval2 = SISL_NULL;
  double *sgt = SISL_NULL;		/* Knot vector of curve put into torus surface      */
  double sy[3];			/* Difference between point and torus center        */
  double tzn;			/* Projection of sy onto torus axis                 */
  double tyy;			/* Square of length of sy                           */
  double tzz;			/* Square of length of sz                           */
  double ty;			/* Component of sy                                  */
  double tz;			/* Component of sz                                  */
  double sder[4];		/* Point on the curve                           */
  double ww;			/* the weight of sder squared if pcurv is rational  */
  double *par = SISL_NULL;
  SISLCurve *tempcurv = SISL_NULL;	/* only used for rational curves              */

  *jstat = 0;
  if (idim != pcurv->idim) goto err104;
  if (ideg != 1001) goto err200;
  
  /* Make local pointers. */

  kn = pcurv->in;
  kk = pcurv->ik;
  kdim = pcurv->idim;
  st = pcurv->et;
  ikind = pcurv->ikind;

  if (ikind == 2 || ikind == 4)
    {
      tempcurv = newCurve (kn, kk, st, pcurv->rcoef, ikind - 1, kdim + 1, 0);
      if (tempcurv == SISL_NULL)
	goto err171;
      tempcurv->cuopen = pcurv->cuopen;
    }
  else
    {
      tempcurv = pcurv;
    }


  /* Test input. */

  if (kdim != 3)
    goto err104;


  /* Make description of knot array for interpolation. */

  s1376 (st, kn, kk, &sgt, &kjkn, &kjkk, &kstat);
  if (kstat < 0)
    goto error;


  /* Make parameter values and derivative indicators. */

  s1890 (sgt, kjkk, kjkn, &par, &der, &kstat);
  if (kstat < 0)
    goto error;


  /* Allocate array for values of curve put into torus equation. */

  sval1 = newarray (kjkn, DOUBLE);
  if (sval1 == SISL_NULL)
    goto err101;


  /* Calculate values to be interpolated. */

  for (ki = 0; ki < kjkn; ki++)
    {
      /*  Calculate values on 3-D curve. */

      s1221 (tempcurv, 0, par[ki], &kleft, sder, &kstat);
      if (kstat < 0)
	goto error;

      /*
       *   The calculation of a point on the torus surface can be done in the
       *   following way.
       *
       *      y = p - scentr
       *      z = y - (y saxis) saxis
       *
       *   The equation of the torus can be written
       *
       *                          2    2
       *      (y - R z/sqrt(z z) )  - r = 0
       *
       *
       *   or by eliminating the square root:
       *
       *          f =
       *
       *          2           2  2      2       2  2 2
       *      (yy)  + 2 (yy)(R -r ) - 4R zz + (R -r )  = 0
       *
       *
       *
       *       or in 4-D homogeneous coordinates:
       *
       *
       *  f =
       *
       *      2      2      2  2      2 2      4  2  2 2
       *  (yy)  + 2 w (yy)(R -r ) - 4w R zz + w (R -r )  = 0
       *
       *      where y = T - w*scentr,  p=T/w
       *
       *   We thus need to calculate yy and zz:
       */

      if (ikind == 2 || ikind == 4)
	{
	  for (kj = 0; kj < 3; kj++)
	    sy[kj] = sder[kj] - sder[3] * scentr[kj];
	  ww = sder[3] * sder[3];
	}
      else
	{
	  for (kj = 0; kj < 3; kj++)
	    sy[kj] = sder[kj] - scentr[kj];
	  ww = (double) 1.0;
	}

      tzn = s6scpr (sy, saxis, 3);

      tyy = (double) 0.0;
      tzz = (double) 0.0;


      /*  Make z and necessary derivatives of z */

      for (kj = 0; kj < 3; kj++)
	{
	  ty = sy[kj];
	  tz = ty - tzn * saxis[kj];
	  tyy += ty * ty;
	  tzz += tz * tz;
	}

/*
 *                                      2            2   2
 * Now tyy = yy and tzz = zz, tbigr2 = R ,tdiffr2 = R - r   
 * --------------------------------------------------------
 */

      sval1[ki] = tyy * tyy + ((double) 2.0 * ww * tyy
			       + ww * ww * tdiffr2) * tdiffr2
	- (double) 4.0 *ww * tbigr2 * tzz;
    }

  cuopen = TRUE;

  s1891 (par, sval1, kone, kjkn, kone, der, cuopen, sgt, &sval2, &kjkn, kjkk,
	 kzero, kzero, &kstat);
  if (kstat < 0)
    goto error;

  *rcurv = SISL_NULL;
  *rcurv = newCurve (kjkn, kjkk, sgt, sval2, 1, 1, 1);
  if (*rcurv == SISL_NULL)
    goto err171;
  (*rcurv)->cuopen = pcurv->cuopen;

 
  /* Ok ! */

  goto out;


  /* Error in lower level function */

error:
  *jstat = kstat;
  s6err ("s1377", *jstat, kpos);
  goto out;

  /* Error in space allocation */

err101:
  *jstat = -101;
  s6err ("s1377", *jstat, kpos);
  goto out;

  /* Dimension not equal to 3 or conflicting dimensions */

err104:
  *jstat = -104;
  s6err ("s1377", *jstat, kpos);
  goto out;

  /* Could not create curve. */

err171:
  *jstat = -171;
  s6err ("s1377", *jstat, kpos);
  goto out;

  /* Wrong implicit type (ideg). */

err200:
  *jstat = -200;
  s6err ("s1377", *jstat, kpos);
  goto out;

out:

  /* Release allocated arrays */

  if (sgt != SISL_NULL)
    freearray (sgt);
  if (par != SISL_NULL)
    freearray(par);
  if (der != SISL_NULL)
    freearray(der);
  if (sval1 != SISL_NULL)
    freearray (sval1);
  if (sval2 != SISL_NULL)
    freearray (sval2);
  if ((ikind == 2 || ikind == 4) && (tempcurv != SISL_NULL))
    freeCurve (tempcurv);

  return;
}
