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
 * $Id: s1733.c,v 1.2 1994-10-19 15:30:55 pfu Exp $
 *
 */


#define S1733

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1733(SISLSurf *ps,int icont1,int icont2,double *cstart1,double *cend1,
	   double *cstart2,double *cend2,double *gcoef,int *jstat)
#else
void s1733(ps,icont1,icont2,cstart1,cend1,cstart2,cend2,gcoef,jstat)
     SISLSurf   *ps;
     int    icont1;
     int    icont2;
     double *cstart1;
     double *cend1;
     double *cstart2;
     double *cend2;
     double *gcoef;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To pick out the next Bezier patch of a B-spline surface.
*              This function requere a B_spline surface that is the
*              result of s1731. This rutine do not check that the
*              surface is correct.
*
*
*
* INPUT      : ps         - B-spline surface to convert.
*              icont1     - The horisontal number of the Bezier patch to pick.
*              icont2     - The vertical number of the Bezier patch to pick.
*
*
*
* OUTPUT     : cstart1    - The horisontal start parameter value to
*                           the Bezier patch.
*              cend1      - The horisontal end parameter value to
*                           the Bezier patch.
*              cstart2    - The vertical start parameter value to
*                           the Bezier patch.
*              cend2      - The vertical end parameter value to
*                           the Bezier patch.
*              gcoef      - The vertices to the Bezier curve.
*              jstat      - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Arne Laksaa, SI, 88-11.
* REVISED BY : Johannes Kaasa, SI, May 1992 (Introduced NURBS).
* REVISED BY : Christophe Birkeland, SI, July 92
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994. Added lower
*              bound checks on 'icont1' and 'icont2' and corrected upper
*              bound checks from "<=" to "<" - caused memory problems.
*
**********************************************************************/
{
  int kpos=0;       /* Position of error.                  */
  int ki;           /* Control variable in loop.           */
  int kfi1,kla1;    /* Index to the first and last element
		       in the knot-vector with the start
		       value to the Bezier curve.          */
  int kfi2,kla2;    /* Index to the first and last element
		       in the knot-vector with the start
		       value to the Bezier curve.          */
  int kdim;         /* Dimension of the geometry.          */
  double *scoef;    /* Vertices.                           */

  *jstat = 0;

  /* Check if the surface is rational. */

  if (ps->ikind == 2 || ps->ikind == 4)
  {
     kdim = ps->idim + 1;
     scoef = ps->rcoef;
  }
  else
  {
     kdim = ps->idim;
     scoef = ps->ecoef;
  }

  /* Check that we have a Bezier patch to treat. */

  if ( icont1 >= 0  &&  icont2 >= 0  &&
       icont1 < ps->in1/ps->ik1  &&  icont2 < ps->in2/ps->ik2 )
    {
      /* Finding the first and last element in ps->et1 with the
	 start1 value. */

      kfi1 = icont1*ps->ik1;
      kla1 = kfi1 + ps->ik1;

      /* Updating the start1 and the end1 parameter value
	 to the patch. */

      *cstart1 = ps->et1[kfi1];
      *cend1 = ps->et1[kla1+1];

      /* Finding the first and last element in ps->et2 with the
	 start2 value. */

      kfi2 = icont2*ps->ik2;
      kla2 = kfi2 + ps->ik2;

      /* Updating the start2 and the end2 parameter value
	 to the patch. */

      *cstart2 = ps->et2[kfi2];
      *cend2 = ps->et2[kla2+1];

      /* Updating the vertices to the Bezier curve. */

      for (ki=0;ki < ps->ik2;ki++)
	memcopy(gcoef+ki*kdim*ps->ik1,
		&scoef[((kfi2+ki)*ps->in1 + kfi1)*kdim],
		kdim*ps->ik1,double);
    }
  else
    {
      /* Error, no patch to return. */

      *jstat = -152;
      s6err("s1733",*jstat,kpos);
    }
}
