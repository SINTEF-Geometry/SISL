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
 * $Id: s1223.c,v 1.5 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1223

#include "sislP.h"
 


#if defined(SISLNEEDPROTOTYPES)
void 
   s1223(double et1[], double et2[], int ik1, int ik2, int in1, int in2, 
	 int ibase1, int ibase2, double par[], int ider1, int ider2, 
	 double ebder[], int *jstat)
#else
void s1223(et1, et2, ik1, ik2, in1, in2, ibase1, ibase2, par, ider1, ider2, 
	   ebder, jstat)
     double et1[];
     double et2[];
     int    ik1;
     int    ik2;
     int    in1;
     int    in2;
     int    ibase1;
     int    ibase2;
     double par[];
     int    ider1;
     int    ider2;
     double ebder[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and derivatives of the tensor product
*              B-spline base function starting at (et1[ibase1], et2[ibase2]), 
*              associated with the knot vectors et1 and et2, at the point par.
*
*
*
* INPUT      : et1    - Double array of dimension [in1+ik1] containing
*                       the knot vector in the first parameter direction.
*              et2    - Double array of dimension [in2+ik2] containing
*                       the knot vector in the second parameter direction.
*              ik1    - The polynomial order of the B-splines associated
*                       with et1.
*              ik2    - The polynomial order of the B-splines associated
*                       with et2.
*              in1    - The dimension of the spline space associated with
*                       the knot vector et1.
*              in2    - The dimension of the spline space associated with
*                       the knot vector et2.
*              ibase1 - The B-spline base function to evaluate (starting at
*                       et1[ibase1]), in the first parameter direction.
*              ibase2 - The B-spline base function to evaluate (starting at
*                       et2[ibase2]), in the second parameter direction.
*              par    - The point at which the tensor product B-spline value 
*                       and derivatives are to be computed.
*              ider1  - The number of derivatives to be computed in the first 
*                       parameter direction.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              ider2  - The number of derivatives to be computed in the second 
*                       parameter direction.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*
*                
*
* OUTPUT     : ebder  - Double array containing value of the tensor product
*                       B-spline and its derivatives at the point par. 
*                       These numbers are stored as rows from the following
*                       matrix:
*                       /                                                 \
*                       | B1(u)B2(v)  DB1(u)B2(v)  ... D(ider1)B1(u)B2(v) |
*                       | B1(u)DB2(v) DB1(u)DB2(v) ...                    |
*                                                 .
*                                                 .
*                       | B1(u)D(ider2)B2(v) ...                          |
*                       \                                                 /
*                       The array has dimension (ider1 + 1)*(ider2 + 1).
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : 
*
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SINTEF, October 1995.
*
*********************************************************************
*/                                     
{
   int kstat=0;          /* Local status variable.                        */
   int kpos=0;           /* The position of the error.                    */
   int ki, kj, kl;       /* Index in for loop.                            */
   int kder1, kder2;     /* Number of derivatives.                        */
   int knumb1;           /* Necessary size of sder1.                      */
   int knumb2;           /* Necessary size of sder2.                      */
   double sdum1[100];    /* Fixed utility array.                          */
   double sdum2[100];    /* Fixed utility array.                          */
   double *sder1 = SISL_NULL; /* Evaluation in the first parameter direction.  */
   double *sder2 = SISL_NULL; /* Evaluation in the second parameter direction. */
   

   /* Establish the number of derivatives. */
   
   kder1 = min(ik1 - 1, ider1);
   kder2 = min(ik2 - 1, ider2);
   
   /* If necessary allocate space for the evaluation in each 
      parameter direction.                                   */
   
   knumb1 = ik1*(kder1 + 1);
   knumb2 = ik2*(kder2 + 1);
   
   if (knumb1 > 100)
   {
      if ((sder1 = newarray(knumb1, double)) == SISL_NULL)
	 goto err101;
   }
   else
      sder1 = &sdum1[0];
   
   if (knumb2 > 100)
   {
      if ((sder2 = newarray(knumb2, double)) == SISL_NULL)
	 goto err101;
   }
   else
      sder2 = &sdum2[0];
   
   /* Compute the value and derivatives in the first parameter direction. */
   
   s1222(et1, ik1, in1, ibase1, par[0], kder1, sder1, &kstat);
   if (kstat < 0) goto error;

   /* Compute the value and derivatives in the second parameter direction. */
   
   s1222(et2, ik2, in2, ibase2, par[1], kder2, sder2, &kstat);
   if (kstat < 0) goto error;
   
   /* Multiply together. */
   
   for (ki = 0, kl = 0; ki < (kder2 + 1); ki++)
   {
      for (kj = 0; kj < (kder1 + 1); kj++, kl++)
      {
	 ebder[kl] = sder1[kj]*sder2[ki];
      }
      for ( ; kj < (ider1 + 1); kj++, kl++)
      {
	 ebder[kl] = 0.;
      }
   }
      
   for ( ; ki < (ider2 + 1); ki++)
   {
      for (kj = 0; kj < (ider1 + 1); kj++, kl++)
      {
	 ebder[kl] = 0.;
      }
   }
   
   /* Successful computations.  */
   
   goto out;
   
   /* Not enough memory. */
   
   err101: 
      *jstat = -101;
   s6err("s1223",*jstat,kpos);
   goto out;
   
   /* Error in lower level routine.  */
   
   error:  
      *jstat = kstat;
   s6err("s1223",*jstat,kpos);
   goto out;
   
   out: 
      if (knumb1 > 100 && sder1 != SISL_NULL) freearray(sder1);
      if (knumb2 > 100 && sder2 != SISL_NULL) freearray(sder2);
      
      return;
}
