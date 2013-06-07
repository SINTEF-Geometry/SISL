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
 * $Id: s1349.c,v 1.3 2001-03-19 15:58:46 afr Exp $
 *
 */


#define S1349

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1349(int inbcrv,SISLCurve *vpcrv[],int *jstat)
#else
void s1349(inbcrv,vpcrv,jstat)
     int   inbcrv;
     SISLCurve *vpcrv[];
     int   *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To convert a B-spline curve with the first ik knots not
*              equal to et(ik) and the last ik knots not equal to
*              et(in+1) to a representation with the first ik knots
*              equal to et(ik), and the last ik knots equal to et(in+1).
*
* INPUT      : inbcrv - No. of curves in the curve-set.
*              vpcrv  - Array (length inbcrv) of pointers to curves
*                       to the curves in the curve-set.
*
* OUTPUT     : jstat  - status messages:
*                        > 0      : warning
*                        = 0      : ok
*                        < 0      : error
*              vpcrv  - Array (length inbcrv) to curves in the curve-set..
*                       (Changed description.)
*
* METHOD     : For each curve, the multiplicity of the end vertices is
*              counted, and those curves not having the right multiplicity
*              are converted to a representation with the right multiplicity.
*-
* CALLS      : s1712,s6err
*
* WRITTEN BY : A. M. Ytrehus SI, Oslo, Norway.  Sep.  1988
* Revised by : Tor Dokken, SI, Oslo, Norway. 26. Feb. 1989
* Revised by : Paal Fugelli, SINTEF, Oslo 02/08-1994. Fixed memory leak.
*
*********************************************************************
*/
{
  SISLCurve **wp = SISL_NULL;      /* Local pointer to current curve.           */
  SISLCurve *qc2 = SISL_NULL;      /* Pointer to new curve-object.              */
  int kvert;              /* No. of vertices in current curve.           */
  int kord;               /* Order of current curve.                     */
  double *sknot = SISL_NULL;   /* Pointer to knot-vector of current curve.    */
  int kk,kr;              /* Loop controllers.                           */
  double *sp1,*sp2;       /* Pointers to sknot.                          */
  int kmul1,kmul2;        /*                                             */
  double tval1,tval2;     /*                                             */
  int kstat = 0;          /* Status variable.                            */
  int kpos = 0;           /* Position of error.                          */

  wp = vpcrv;
  for (kk=0; kk<inbcrv; kk++)
    {
      /* Make local pointers to description of curve */
      sknot = (*wp) -> et;
      kvert = (*wp) -> in;
      kord = (*wp) -> ik;

      /* Count multiplicity of start-knot. */
      kmul1 = 0;
      sp1 = sknot + kord - 1;
      tval1 = *sp1;
      for (kr=0; kr<kord; kr++)
	{
          if (*sp1 == tval1) kmul1++;
          sp1--;
	}

      /* Count multiplicity of end-knot. */
      kmul2 = 0;
      sp2 = sknot + kvert;
      tval2 = *sp2;
      for (kr=0; kr<kord; kr++)
	{
          if (*sp2 == tval2) kmul2++;
          sp2++;
	}

      /* If the multiplicity of both end-knots equals kord, the curve is ok. */
      if (kmul1 != kord || kmul2 != kord)
	{
          /* Both ends do not have multiplicity of order kord.
             Create a new curve-object.                        */
          s1712((*wp),tval1,tval2,&qc2,&kstat);
          if (kstat<0) goto error;
	  if ((*wp)) freeCurve(*wp);  /* PFU 02/08-1994 */
          *wp = qc2;
          qc2 = SISL_NULL;
	}
      wp++;
    }

  *jstat = 0;
  goto out;

/* Error in lower level routine. */
error: *jstat = kstat;
       s6err("s1349",*jstat,kpos);
       goto out;
out:
return;
}
