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
 * $Id: s1959.c,v 1.3 2001-03-19 15:58:58 afr Exp $
 *
 */
#define S1959

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1959(SISLPoint *ppoint, SISLCurve *pcurve, double *gpos, int *jstat)
#else
void s1959(ppoint,pcurve,gpos,jstat)
     SISLPoint        *ppoint;
     SISLCurve         *pcurve;
     double       *gpos;
     int          *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Estimate parameter of guess-point (used by closest point
*              calculation).
*
*
* INPUT      : ppoint   - Pointer to the point.
*              pcurve    - Pointer to the curve.
*
* OUTPUT     : gpos    - Parameter of the found guess-point.
*              jstat   - status messages
*                                = 0   : Guess-point found.
*                                = 1   : Closest point as guess-point.
*                                < 0   : error.
*
*
* METHOD     : Quadrant analysis and Schoenberger equation.
*
*
* REFERENCES :
*
*
* WRITTEN BY : Michael Floater, SI, October 1991
*                Based on Per Evensen's surface version s1960.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Dec. 1994. Fixed memory
*              problem for 2D curves.
*
*********************************************************************
*/
{
  int kstat = 0;            /* Local status variable.                     */
  int kpos = 0;             /* Position of error.                         */
  int i;                    /* Running indexes                            */
  int iind;                 /* Index variable                             */
  int kk;                /* The polynomial order of the curve (pcurve) */
  int nk;                /* The number of vertices of the curve (pcurve)  */
  int dim;                  /* Dimension of space the curve lies in      */
  double *et;              /* Knots for the curve (pcurve)               */
  double *lpoint;           /* Coefficients of the point (ppoint)         */
  double *scoef;            /* Coefficients of the curve (pcurve)         */
  double tdist;             /* Distance variable                          */
  double tmin;              /* Minimum distance variable                  */
  double vec1[3],vec2[3];   /* Vectors defining the quadrants surrounding
                               a vertex                                  */
  double vecp[3];           /* Relative point vector                      */
  double lqua1=DZERO;
  double lqua2=DZERO;
                            /* Length of vec1 and vec2                   */
  double lprj1=DZERO;
  double lprj2=DZERO;
                            /* Length of projection of vecp on
                               vec1 and vec2                             */
  double svals,svale;       /* Local start and end parameter values      */
 /* --------------------------------------------------------------------- */

  /* Test input.  */
  if (ppoint->idim != pcurve->idim || pcurve->idim <= 1) goto err106;

  /* Initialize local variables.  */
  kk  = pcurve->ik;
  nk  = pcurve->in;
  et = pcurve->et;
  scoef = pcurve->ecoef;
  dim = pcurve->idim;
  lpoint = ppoint->ecoef;

  /* Find vertex closest to point.  */
  tdist=s6dist(scoef,lpoint,dim);
  tmin=tdist;
  iind = 0;
  for (i=0; i<nk; i++)
  {
        tdist=s6dist(scoef,lpoint,dim);
        if (tdist<tmin)
        {
           tmin=tdist;
           iind = i;
        }
        scoef+=dim;  /* PFU 02/12-94 Changed from 'scoef+=3' */
  }

  /*
  * Generate the "vectors of the quadrant".



           vec2              vec 1
             <-------X------->
closest            /     x
     vertex (iind)        \ point

  */
  scoef = pcurve->ecoef;

  if (iind < (nk-1))
     s6diff(&(scoef[(iind+1)*dim]),
            &(scoef[iind*dim]),dim,vec1);
  if (iind > 0)
     s6diff(&(scoef[(iind-1)*dim]),
            &(scoef[iind*dim]),dim,vec2);

  /* Generate the point - closest vertex vector. */
  s6diff(lpoint,&(scoef[iind*dim]),dim,vecp);

  /* Calculate the length of the quadrant vectors. */
  if (iind < (nk-1)) lqua1 = s6length(vec1,dim,&kstat);
  if (iind > 0) lqua2 = s6length(vec2,dim,&kstat);

  /* Calculate the length of the projection of 'vecp' on the quadrant
     vectors. */
  if (iind < (nk-1)) lprj1 = s6lprj(vecp,vec1,dim);
  if (iind > 0) lprj2 = s6lprj(vecp,vec2,dim);

  /* Decide in which quadrant the point lies. */
  if (iind == 0)
  {
     /* Point lies in 1. quadrant */

     /* Calculate knot values of vertex */
     svals = s6schoen(et,kk,iind);
     svale = s6schoen(et,kk,iind+1);

     /* Calculate estimated parameter values of point */
     if (lqua1 != DZERO) (*gpos) = svals + (lprj1/lqua1)*(svale-svals);
     else (*gpos) = svals;
  }
  else if (iind == (nk-1))
  {
     /* Point lies in 2. quadrant */

     /* Calculate knot values of vertices */
     svals = s6schoen(et,kk,iind-1);
     svale = s6schoen(et,kk,iind);

     /* Calculate estimated parameter values of point */
     if (lqua2 != DZERO) (*gpos) = svals + ((lqua2-lprj2)/lqua2)*(svale-svals);
     else (*gpos) = svals;
  }
  else if (iind > 0 && iind < (nk-1))
  {
     /* Evaluate 1. and 2. quadrant */

      if (lprj1 > lprj2)
      {
           /* Point lies in 1. quadrant */

           /* Calculate knot values of vertices */
           svals = s6schoen(et,kk,iind);
           svale = s6schoen(et,kk,iind+1);

           /* Calculate estimated parameter values of point */
           if (lqua1 != DZERO) (*gpos) = svals + (lprj1/lqua1)*(svale-svals);
           else (*gpos) = svals;
      }
      else if (lprj2 > lprj1)
      {
           /* Point lies in 2. quadrant */

           /* Calculate knot values of vertices */
           svals = s6schoen(et,kk,iind-1);
           svale = s6schoen(et,kk,iind);

           /* Calculate estimated parameter values of point */
           if (lqua2 != DZERO) (*gpos) = svals + ((lqua2-lprj2)/lqua2)*(svale-svals);
           else (*gpos) = svals;
      }
      else
      {
           /* lprj1 and lprj2 are both zero or equal. */
           /* Choose the control point itself. */

           goto usvert;
      }
  }
  else
  {
     /* Error */
     goto usvert;
  }

  /* Check that values are within parameter plane of curve. */
  if ((*gpos)<et[kk-1]) (*gpos)=et[kk-1];
  else if ((*gpos)>et[nk]) (*gpos)=et[nk];
  *jstat = 0;

  /* Calculation completed.  */
  goto out;

  /* No intermediate parameter value found,
     use parameter value of closest vertex. */
  usvert: *jstat = 1;

  /* Calculate knot value of closest vertex */
  (*gpos) = s6schoen(et,kk,iind);

  /* Check that values are within parameter plane of curve. */
  if ((*gpos)<et[kk-1]) (*gpos)=et[kk-1];
  else if ((*gpos)>et[nk]) (*gpos)=et[nk];

  goto out;

 /* --------------------------------------------------------------------- */
  /* Error in input. Dimension not equal to 1 */
 err106: *jstat = -106;
  s6err("s1959",*jstat,kpos);
  goto out;

 out:
    return;
}
