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

#define S1522

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1522(double normal[], double centre[], double ellipaxis[],
           double ratio, int dim, SISLCurve** ellipse, int *jstat)
#else
void s1522(normal, centre, ellipaxis, ratio, dim, ellipse, jstat)
     double normal[];
     double centre[];
     double ellipaxis[];
     double ratio;
     int dim;
     SISLCurve** ellipse;
     int *jstat;
#endif
/*
***************************************************************************
*
*  PURPOSE      : Convert a 2D or 3D Ellipse to non-uniform rational B-spline
*                 (nurbs) format.
*
*  INPUT        :
*        normal - 3D normal to ellipse plane (not necessarily normalized).  Used
*                 if dim = 3.
*        centre - centre of ellipse (2D if dim = 2 and 3D if dim = 3).
*     ellipaxis - One axis of the ellipse -- will be used as starting point
*                 for the ellipse curve (2D if dim = 2 and 3D if dim = 3).
*         ratio - The ratio between the length of the given ellipaxis and the
*                 length of the other axis, i.e. |ellipaxis| / |other axis| (a
*                 compact representation format).
*           dim - dimension of the space in which the elliptic nurbs curve
*                 lies (2 or 3).
*
*  OUTPUT       :
*   ellipse - ellipse in nurbs format (2D if dim = 2 and 3D if dim = 3).
*   jstat   - status messages
*                         > 0      : warning
*                         = 0      : ok
*                         < 0      : error
*
*  METHOD       : The length of the parameter interval is approximately
*                 equal to the length of the ellipse so that the
*                 parametrization is approximately arc length.
*
*  CALLS        : newCurve().
*
*
*  REFERENCES   :
*
*
* WRITTEN BY    :  Michael Floater, based on earlier version by
*                  Paul Fugelli.      SINTEF SI (Norway), August 1994.
*
****************************************************************************/
{
  int kstat;       /* Status variable.                    */
  int kpos=0;      /* Position of error.                  */
  double cross[3]; /* Cross vector (the other axis). */
  int i;           /* For-loop indices.   */
  double et[12];   /* Knot vector.        */
  double coef[36]; /* Vertices.           */
  double elliplen; /* Length of ellipaxis. */
  double otherlen; /* Length of the other axis. */
  double factor;   /* Length multiplication factor. */
  double tworoot = sqrt ((double) 2.0);  /* Square root of 2. */
  double weight  = (double) 1.0 / tworoot;  /* One over square root of 2. */

  *jstat = 0;

  /* Check input. */
  if ( ratio == 0.0 )  ratio = 1.0;

  if ( dim != 2  &&  dim != 3 )  goto err105;


  /* Find the length of the given axis. */

  elliplen = s6length(ellipaxis, dim, &kstat);
  if ( kstat < 0  ||  elliplen == 0.0 )  goto error;

  /* Find the other axis. */
  if ( dim == 2 )
  {
    cross[0] = -ellipaxis[1] / ratio;
    cross[1] = ellipaxis[0] / ratio;

    /* Find the length of the other axis. */
    otherlen = elliplen / ratio;
  }
  else
  {
    cross[0] = normal[1]*ellipaxis[2] - normal[2]*ellipaxis[1];
    cross[1] = normal[2]*ellipaxis[0] - normal[0]*ellipaxis[2];
    cross[2] = normal[0]*ellipaxis[1] - normal[1]*ellipaxis[0];

    /* Find the length of the other axis. */
    otherlen = s6length(cross, 3, &kstat);
    if ( kstat < 0  ||  otherlen == 0.0 )  goto error;

    /* Normalise the length of the other axis. */
    for ( i=0;  i < 3;  i++ )  cross[i] /= otherlen;

    /* Find the length of the other axis. */
    otherlen = elliplen / ratio;

    /* Adjust the length of the other axis. */
    for ( i=0;  i < 3;  i++ )  cross[i] *= otherlen;
  }


  /* Build the nurbs representation. */

  factor = PI*tworoot * sqrt( elliplen*elliplen + otherlen*otherlen );
  et[0] = (DOUBLE) 0.0;  /* The knot vector. */
  for ( i=1;  i < 3;  i++ ) {
    et[i]     = (DOUBLE) 0.0;
    et[2 + i] = factor * 0.25;
    et[4 + i] = factor * 0.5;
    et[6 + i] = factor * 0.75;
    et[8 + i] = factor;
  }
  et[11] = factor;

  if ( dim == 2 )
  {
    for ( i=0;  i < 2;  i++ )
    {
      coef[     i] = centre[i] + ellipaxis[i];
      coef[3 +  i] = weight*(centre[i] + ellipaxis[i] + cross[i]);
      coef[6 +  i] = centre[i] + cross[i];
      coef[9 + i] = weight*(centre[i] - ellipaxis[i] + cross[i]);
      coef[12 + i] = centre[i] - ellipaxis[i];
      coef[15 + i] = weight*(centre[i] - ellipaxis[i] - cross[i]);
      coef[18 + i] = centre[i] - cross[i];
      coef[21 + i] = weight*(centre[i] + ellipaxis[i] - cross[i]);
      coef[24 + i] = centre[i] + ellipaxis[i];
    }

    coef[2] = 1.0;  /* The rational weights. */
    coef[5] = weight;
    coef[8] = 1.0;
    coef[11] = weight;
    coef[14] = 1.0;
    coef[17] = weight;
    coef[20] = 1.0;
    coef[23] = weight;
    coef[26] = 1.0;
  }
  else
  {
    for ( i=0;  i < 3;  i++ )
    {
      coef[     i] = centre[i] + ellipaxis[i];
      coef[4 +  i] = weight*(centre[i] + ellipaxis[i] + cross[i]);
      coef[8 +  i] = centre[i] + cross[i];
      coef[12 + i] = weight*(centre[i] - ellipaxis[i] + cross[i]);
      coef[16 + i] = centre[i] - ellipaxis[i];
      coef[20 + i] = weight*(centre[i] - ellipaxis[i] - cross[i]);
      coef[24 + i] = centre[i] - cross[i];
      coef[28 + i] = weight*(centre[i] + ellipaxis[i] - cross[i]);
      coef[32 + i] = centre[i] + ellipaxis[i];
    }

    coef[3] = 1.0;  /* The rational weights. */
    coef[7] = weight;
    coef[11] = 1.0;
    coef[15] = weight;
    coef[19] = 1.0;
    coef[23] = weight;
    coef[27] = 1.0;
    coef[31] = weight;
    coef[35] = 1.0;
  }

  (*ellipse) = newCurve(9, 3, et, coef, 2, dim, 1);
  if ( (*ellipse) == SISL_NULL )  goto err101;

  (*ellipse)->cuopen = 0;

  goto out;


  /* Error in curve allocation.  */

err101:
  *jstat = -101;
  s6err("s1522", *jstat, kpos);
  goto out;


  /* Error in input. Dimension not equal to 2 or 3. */

err105:
  *jstat = -105;
  s6err("s1522", *jstat, kpos);
  goto out;


  /* Error in lower level routine. */

error:
  *jstat = kstat;
  s6err("s1522", *jstat, kpos);
  goto out;


out:
  return;
}
