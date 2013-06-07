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
 * $Id: s1391.c,v 1.3 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1391

#include "sislP.h"

typedef void (*fshapeProc)(
#if defined(SISLNEEDPROTOTYPES)
                           double [],
                           double [],
                           int,
                           int,
                           int *
#endif
);


#if defined(SISLNEEDPROTOTYPES)
void
s1391(SISLCurve **pc,SISLSurf ***ws,int icurv,int nder[],int *jstat)
#else
void s1391(pc,ws,icurv,nder,jstat)
     SISLCurve **pc;
     SISLSurf  ***ws;
     int   icurv;
     int   nder[];
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* Purpose   : To create a first derivative continuous blend over a
*             3-, 4-, 5- and 6-sided region in space. The boundary of the
*             region are B-spline curves and the cross boundary derivatives
*             are given as B-spline curves. This function automatically
*             preprosesses the input cross tangent curves in order to
*             make them suitable for the blending. Thus, the cross tangent
*             curves should be taken as the cross tangents of the
*             surrounding surface. It is not necessary and not advisable
*             to match tangents etc. in the corners.
*
*
* Input     : pc        - Pointers to boundary curves. All curves must
*                         have same parameter direction around the patch,
*                         either clockwise or counter-clockwise.
*                         pc1[i],i=0,...nder[0]-1 are pointers to position
*                         and cross-derivatives along first edge.
*                         pc1[i],i=nder[0],...nder[1]-1 are pointers
*                         to position and cross-derivatives along second edge.
*                         .
*                         .
*                         pc1[i],i=nder[0]+...+nder[icurve-2],...
*                         nder[icurve-1]-1 are pointers to position and
*                         cross-derivatives along fourth edge.

*             icurv     - (3,5,4 or 6), Number of boundary curves.
*             nder[0:icurv-1] - nder[i] gives number of curves on
*                               edge number i+1. These numbers has to
*                               be equal to 2.
*
* Output    : ws        - ws[0:icurve-1] are pointers to the blending surface
*
*             jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* Use        : SISLCurve *pc1[KCURVE];
*              SISLSurf  **rs1;
*              int jstat,icurve,nder[KSURF];
*
*              icurve = number of boundaries;
*              nder[0] = nder[1] = ... = nder[icurve-1] = 2;
*              pc1[0] = curve one;
*              .
*              .
*
*              s1391(pc1,&rs1,icurve,nder,&jstat);  icurve surfaces is given
*                                                   as output.
*
*
*  NB!         Position and tangent has to be given, that is
*              n1=n2=n3=((n4=n5)=n6)=2.
*
*-
* Calls      : s1720  - Differentiate B-spline curve.
*              sh1263 - Prepare input curves for the blending
*              sh1460 - Perform blending.
*              s6err  - Error messages.
*
* Written by : Morten Daehlen, SI, Aug. 88.
* Rewritten by : Vibeke Skytt, SINTEF SI, 93-11.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994. Fixed memory
*              problems related to 'qc1'.
*
*********************************************************************
*/
{
   int kstat = 0;   /* Local status variable.  */
   int kpos = 0;
   int ki;
   SISLCurve **qc1 = SISL_NULL;  /* Input curves to the
			       preparation of curves function.         */
   SISLCurve **qc2 = SISL_NULL;  /* Output curves from the
			       preparation of curves function.         */
   fshapeProc fshape;
   /* #if defined(SISLNEEDPROTOTYPES)
   void (*fshape)(double*, double*, int, int, int*);
#else
   void (*fshape)();
#endif  */
   fshape = shape;

   if ((qc1 = newarray(3*icurv, SISLCurve*)) == SISL_NULL) goto err101;
   if ((qc2 = newarray(2*icurv, SISLCurve*)) == SISL_NULL) goto err101;
   memzero(qc1, 3*icurv, SISLCurve*);
   memzero(qc2, 2*icurv, SISLCurve*);

   /* Test input.  */

   if (icurv != 3 && icurv != 4 && icurv != 5 && icurv != 6) goto err180;
   for (ki=0; ki<icurv; ki++)
      if (nder[ki] != 2) goto err180;

   /* Copy input curve to the prepare array. Differentiate
      the position curve to get a second derivative curve. */

   for (ki=0; ki<icurv; ki++)
   {
      qc1[3*ki] = pc[2*ki];
      qc1[3*ki+1] = pc[2*ki+1];
      qc1[3*ki+2] = SISL_NULL;

      s1720(qc1[3*ki], 1, qc1+3*ki+2, &kstat);
      if (kstat < 0) goto error;
   }

   /* Prepare boundary conditions for vertex blend. */

   sh1263(qc1, icurv, qc2, &kstat);
   if (kstat < 0) goto error;

   /* Perform blending of 3-, 4-, 5 or 6-sided area.  */

   sh1460(fshape, qc2, icurv, ws, &kstat);
   if (kstat < 0) goto error;

   /* Blending performed. */

   *jstat = 0;
   goto out;

   /* Error in scratch allocation. */

   err101 :
      *jstat = -101;
   s6err("s1391", *jstat, kpos);
   goto out;

   /* Error in input. Wrong number of boundaries or derivative curves. */

   err180 :
      *jstat = -180;
   s6err("s1391", *jstat, kpos);
   goto out;

   /* Error in lower level function. */

   error :
      *jstat = kstat;
   s6err("s1391", *jstat, kpos);
   goto out;

   out :
      for (ki=0; ki<icurv; ki++)
      {
	 if (qc1[3*ki+2] != SISL_NULL) freeCurve(qc1[3*ki+2]);
	 if (qc2[2*ki] != SISL_NULL) freeCurve(qc2[2*ki]);
	 if (qc2[2*ki+1] != SISL_NULL) freeCurve(qc2[2*ki+1]);
      }
      if (qc1 != SISL_NULL) freearray(qc1);
      if (qc2 != SISL_NULL) freearray(qc2);

      return;
}
