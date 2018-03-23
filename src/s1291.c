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


#define S1291

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1291(double cvder[], double sfder[], int idim, double cvder2[], 
      int *jstat)
#else
  void s1291(cvder, sfder, idim, cvder2, jstat)
     double cvder[];
     double sfder[];
     int idim;
     double cvder2[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate position first and second derivative of
*              the projection of a curve in a surface to a corresponding
*              curve
*
*
*
* INPUT      : cvder   - 0-2 order derivatives of the corresponding curve
*              sfder   - 0-2 order derivatives of surface. 
*                        The sequence is position, first derivative in first
*                        parameter direction, first derivative in second
*                        parameter direction, (2,0) derivative, (1,1)
*                        derivative, (0,2) derivative and normal. (21 numbers)
*                        Compatible with output of s1421
*              idim    - Dimension of geometry space, expected to be 2 or 3
*
* OUTPUT     : cvder2  - 0-2 order derivatives of curve in surface
*              jstat  - status messages  
*                         = 0      : ok
*                         < 0      : error
*
* METHOD     : The 0-2 order derivatives of the parameter curve in the
*              surface is computed by projection. Then the derivatives
*              of the surface curve is computed using the chain rule
*
*
* REFERENCES : 
*-
*              
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 02.2018
*
*********************************************************************
*/
{
  int ki;                     /* Counter                            */
  double *f_t = cvder+idim;  /* 1. derivative, curve               */
  double *f_tt = f_t+idim;  /* 2. derivative, curve               */
  double *s_u = sfder+idim;  /* Partial derivative, 1. par. dir surface */
  double *s_v = s_u+idim;   /* Partial derivative, 2. par. dir surface */
  double *s_uu = s_v+idim;  /* 2. order partial derivative, 1. par. dir. */
  double *s_uv = s_uu+idim; /* Mixed 2. order derivative          */
  double *s_vv = s_uv+idim; /* 2. order partial derivative, 2. par. dir. */
  double susu, susv, svsv, suft, svft;  /* Scalar products between 1.
					   order derivatives               */
  double tdiv;                /* Divisor in equation solutions             */
  double c1 = 0.0, c2 = 0.0;  /* 1. order derivative of curve in 
				 parameter domain                   */
  double c3 = 0.0, c4 = 0.0;  /* 1. order derivative of curve in 
				 parameter domain                   */
  double v1[3], v2[3];        /* Help vectors in computations       */
  double t2u, t2v;            /* Scalar products                    */

  /* Check input */
  if (idim != 2 && idim != 3)
    goto err105;

  /* Compute 1. order derivative with respect to the parameter domain curve */
  susu = s6scpr(s_u, s_u, idim);
  susv = s6scpr(s_u, s_v, idim);
  svsv = s6scpr(s_v, s_v, idim);
  suft = s6scpr(s_u, f_t, idim);
  svft = s6scpr(s_v, f_t, idim);
  tdiv = susv*susv - susu*svsv;
  if (fabs(tdiv) < REL_PAR_RES)
    {
      if (fabs(susu) < REL_PAR_RES && fabs(svsv) < REL_PAR_RES)
	{
	  c2 = s6scpr(f_t, f_t, idim)/svsv;
	  if (s6scpr(s_v, f_t, idim) < 0.0)
	    c2 *= -1.0;
	}
      else
	{
	  c1 = s6scpr(f_t, f_t, idim)/susu;
	  if (s6scpr(s_u, f_t, idim) < 0.0)
	    c1 *= -1.0;
	}
    }
  else
    {
      c1 = (svft*susv - suft*svsv)/tdiv;
      c2 = (suft*susv - svft*susu)/tdiv;
    }

  /* Compute 2. order derivative with respect to the parameter domain curve */
  for (ki=0; ki<idim; ++ki)
    {
      v1[ki] = s_uu[ki]*c1*c1 + 2*s_uv[ki]*c1*c2 + s_vv[ki]*c2*c2;
      v2[ki] = f_tt[ki] - v1[ki];
    }
  
  t2u = s6scpr(v2, s_u, idim);
  t2v = s6scpr(v2, s_v, idim);
  if (fabs(tdiv) >= REL_PAR_RES)
    {
      c3 = (t2u*susv - t2v*svsv)/tdiv;
      c4 = (t2u*susv - t2v*susu)/tdiv;
    }

  /* Compute derivatives with respect to the surface curve */
  for (ki=0; ki<idim; ++ki)
    {
      cvder2[ki] = sfder[ki];
      cvder2[idim+ki] = s_u[ki]*c1 + s_v[ki]*c2;
      cvder2[2*idim+ki] = s_uu[ki]*c1*c1 + s_u[ki]*c3 + s_vv[ki]*c2*c2 +
	s_v[ki]*c4 + 2*s_uv[ki]*c1*c2;
    }

  goto out;

 err105:
  *jstat = -105;
  goto out;

 out:
  return;
}
