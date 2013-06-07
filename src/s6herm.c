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
 * $Id: s6herm.c,v 1.3 2005-02-28 09:04:50 afr Exp $
 *
 */


#define S6HERM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s6herm(double *pt,double *uknots,double *vknots,int unum,int vnum,
	 int dim,int uindex,int vindex,double herminfo[],int *jstat)
#else
void
s6herm(pt,uknots,vknots,unum,vnum,dim,uindex,vindex,herminfo,jstat)
     double *pt;
     double *uknots;
     double *vknots;
     int    unum;
     int    vnum;
     int    dim;
     int    uindex;
     int    vindex;
     double herminfo[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Fit p_u, p_v derivatives and cross deriv p_uv to
*              one point in a given rectangular grid of interpolation
*              points.
*              The derivatives are taken from the
*              local biquadratic surface which interpolates the
*              nearest 9 points. This information can be used to
*              fit a cubic Hermite surface.
*
*
*
* INPUT      : pt     - The input grid of points.
*                       This should contain unum*vnum*dim doubles.
*              uknots - The parameter values of the points in u direc.
*              vknots - The parameter values of the points in u direc.
*              unum   - Number of points (and knots) in u direc (>=2).
*              vnum   - Number of points (and knots) in u direc (>=2).
*              dim    - Dimension of points (<= 3).
*              uindex - Index in u of given point (0 to unum-1).
*              vindex - Index in v of given point (0 to unum-1).
*
*
* OUTPUT     : herminfo - The derivs p_u, p_v, p_uv. This
*                         array must have room for 9 doubles
*                         (3 vectors of max dimension 3).
*              jstat  - status messages
*                     = 2      : Surface is degenerate
*                                at the point, normal
*                                has zero length
*                     = 1      : Surface is close to
*                                degenerate at the point
*                                Angle between tangents,
*                                less than angular tolerance
*                     = 0      : ok
*                     < 0      : error
*
* METHOD     : The unique interpolating biquadratic polynomial
*              which fits the nearest 3x3 grid of points is
*              calculated (if unum > 2, vnum > 2).
*              The derivatives are then found
*              from this function.
*              The function is expressed explicitly by
*              Lagrange polynomials.
*              If unum = 2, we fit a 1 x 2 surface.
*              If vnum = 2, we fit a 2 x 1 surface.
*              If unum = vnum = 2, we fit a 1 x 1 surface.
*
* EXAMPLE OF USE:
*
*
* REFERENCES :
*
*
*
* WRITTEN BY : Michael Floater, SI, 3/2/92.
*
*********************************************************************
*/
{
  int kpos = 0;      /* Position of error.          */
  int ki,kj,kk;       /* Loop variable. */
  double Lu[3];         /* Lagrange polynomials in u. */
  double Lv[3];         /* Lagrange polynomials in v. */
  double Ldu[3];        /* 1st derivatives of Lagrange polys in u. */
  double Ldv[3];        /* 1st derivatives of Lagrange polys in v. */
  double temp[3],temp2[3]; /* Temporary storage. */
  double upar,vpar;        /* Parameter values to evaluate at. */
  double diff[3];     /* Temporary storage for Lag polys. */
  double kdiff[3];    /* Temporary storage for Lag polys. */
  int uind,vind;   /* Indices of umin,vmin, in interpolation grid. */
  int i00;                  /* Index in pt array of umin,vmin. */
  int udeg,vdeg; /* Degree of polynomial in each direction (1 or 2). */
  int ind;  /* Temporary index. */

  /* Check input. */

  if(dim < 1 || dim > 3) goto err105;
  if(unum < 2  || vnum < 2) goto err105;
  if(uindex < 0  || uindex > unum) goto err105;
  if(vindex < 0  || vindex > vnum) goto err105;

  /* Set parameter values to evaluate biquadratic at. */

  upar = uknots[uindex];
  vpar = vknots[vindex];

  /* Decide on the degree of the interpolating polynomial. */

  udeg = (unum > 2 ? 2 : 1);
  vdeg = (vnum > 2 ? 2 : 1);

  /* Find bottom left hand corner of grid to interpolate at. */

  if(udeg == 2)
  {
      if(uindex == 0) uind = 0;
      if(uindex > 0) uind = uindex - 1;
      if(uindex == unum-1) uind = uindex - 2;
  }
  else
  {
      uind = 0;
  }

  if(vdeg == 2)
  {
      if(vindex == 0) vind = 0;
      if(vindex > 0) vind = vindex - 1;
      if(vindex == vnum-1) vind = vindex - 2;
  }
  else
  {
      vind = 0;
  }


  /* Calculate Lagrange polynomials in u
     and calculate their 1st derivatives. */

  if(udeg == 2)
  {
      diff[0] = upar-uknots[uind];
      diff[1] = upar-uknots[uind+1];
      diff[2] = upar-uknots[uind+2];

      kdiff[0] = uknots[uind]-uknots[uind+1];
      kdiff[1] = uknots[uind]-uknots[uind+2];
      kdiff[2] = uknots[uind+1]-uknots[uind+2];

      Lu[0] =   diff[1] * diff[2] / (kdiff[0] * kdiff[1]);
      Lu[1] = - diff[0] * diff[2] / (kdiff[0] * kdiff[2]);
      Lu[2] =   diff[0] * diff[1] / (kdiff[1] * kdiff[2]);

      Ldu[0] =   (diff[1] + diff[2]) / (kdiff[0] * kdiff[1]);
      Ldu[1] = - (diff[0] + diff[2]) / (kdiff[0] * kdiff[2]);
      Ldu[2] =   (diff[0] + diff[1]) / (kdiff[1] * kdiff[2]);
  }
  else
  {
      diff[1] = upar-uknots[uind];
      diff[2] = upar-uknots[uind+1];

      kdiff[2] = uknots[uind]-uknots[uind+1];

      Lu[0] =   diff[2]  / kdiff[2];
      Lu[1] = - diff[1]  / kdiff[2];
      Lu[2] =   DZERO;

      Ldu[0] =   (double)1.0  / kdiff[2];
      Ldu[1] = - (double)1.0  / kdiff[2];
      Ldu[2] = DZERO;
  }


  /* Calculate Lagrange polynomials in v
     and calculate their 1st derivatives. */

  if(vdeg == 2)
  {
      diff[0] = vpar-vknots[vind];
      diff[1] = vpar-vknots[vind+1];
      diff[2] = vpar-vknots[vind+2];

      kdiff[0] = vknots[vind]-vknots[vind+1];
      kdiff[1] = vknots[vind]-vknots[vind+2];
      kdiff[2] = vknots[vind+1]-vknots[vind+2];

      Lv[0] =   diff[1] * diff[2] / (kdiff[0] * kdiff[1]);
      Lv[1] = - diff[0] * diff[2] / (kdiff[0] * kdiff[2]);
      Lv[2] =   diff[0] * diff[1] / (kdiff[1] * kdiff[2]);

      Ldv[0] =   (diff[1] + diff[2]) / (kdiff[0] * kdiff[1]);
      Ldv[1] = - (diff[0] + diff[2]) / (kdiff[0] * kdiff[2]);
      Ldv[2] =   (diff[0] + diff[1]) / (kdiff[1] * kdiff[2]);
  }
  else
  {
      diff[1] = vpar-vknots[vind];
      diff[2] = vpar-vknots[vind+1];

      kdiff[2] = vknots[vind]-vknots[vind+1];

      Lv[0] =   diff[2]  / kdiff[2];
      Lv[1] = - diff[1]  / kdiff[2];
      Lv[2] =   DZERO;

      Ldv[0] =   (double)1.0  / kdiff[2];
      Ldv[1] = - (double)1.0  / kdiff[2];
      Ldv[2] = DZERO;
  }

  /* Calculate derivative of biquadratic in u. */

  i00 = (vind*unum+uind)*dim;

  for(ki=0; ki<dim; ki++)
  {
      ind = i00+ki;

      for(kj=0; kj<=vdeg; kj++)
      {
	  temp[kj] = 0.0;
	  temp2[kj] = 0.0;

          for(kk=0; kk<=udeg; kk++,ind+=dim)
          {
	      temp[kj] += Lu[kk] * pt[ind];
	      temp2[kj] += Ldu[kk] * pt[ind];
	  }

	  ind += (unum - udeg - 1) * dim;
      }

      herminfo[ki] = 0.0;
      herminfo[dim+ki] = 0.0;
      herminfo[dim+dim+ki] = 0.0;

      for(kj=0; kj<=vdeg; kj++)
      {
	  /* Calculate u derivative. */

	  herminfo[ki] += Lv[kj] * temp2[kj];

	  /* Calculate v derivative. */

	  herminfo[dim+ki] += Ldv[kj] * temp[kj];

	  /* Calculate uv derivative. */

	  herminfo[dim+dim+ki] += Ldv[kj] * temp2[kj];
      }
  }



  /* Derivatives calculated. */

  *jstat = 0;
  goto out;

  /* Error in input. */

  err105: *jstat = -105;
  s6err("s6herm",*jstat,kpos);
  goto out;

  out: return;
}
