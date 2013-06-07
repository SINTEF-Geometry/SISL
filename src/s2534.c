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
 * $Id: s2534.c,v 1.2 2001-03-19 15:58:59 afr Exp $
 *
 */


#define S2534

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s2534(SISLSurf *surf, 
	  int u_multinc, 
	  int v_multinc, 
	  int newik1, 
	  int newik2,
	  void evalp(SISLSurf *surf, int ider, int iside1, int iside2, 
		     double parvalue[], int *leftknot1, int *leftknot2, 
		     double *result, int *istat),
	  int eval_dim, 
	  SISLSurf **rsurf, 
	  int *stat)
#else
   void 
      s2534(surf, u_multinc, v_multinc, newik1, newik2, evalp, eval_dim, 
	     rsurf, stat)
      SISLSurf *surf;
      int u_multinc;
      int v_multinc;
      int newik1; 
      int newik2;
      void evalp();
      int eval_dim;
      SISLSurf **rsurf;
      int *stat;
#endif
/*
*********************************************************************
*
* PURPOSE : To derive a properity surface from an original surface.
*           The new spline space is based on the original, but new orders
*           and knot multiplicities are given. The property is evaluated with
*           the evaluator evalp.
*
*           We assume that the input knot vectors are k-regular, and that the
*           knot multiplicity is prechecked to avoid interior knot multiplicity
*           equal to or larger than the order.
*
*
*
* INPUT   : surf      - The original k-regular surface.
*           u_multinc - The multiplicity increment in the first direction.
*                       In addition the multiplicity is increased by the order
*                       increase.
*           v_multinc - The multiplicity increment in the second direction.
*                       In addition the multiplicity is increased by the order
*                       increase.
*           newik1    - The new order in the first direction.
*           newik2    - The new order in the second direction.
*           evalp     - The generic property evaluator.
*           eval_dim  - Dimension of the result from the evaluator.
*
*
*
* OUTPUT   : rsurf    - The resulting surface
*            stat     - Status messages
*                       > 0      : Warning
*                       = 2      : Degenerated surface
*                       = 0      : Ok
*                       < 0      : Error
*
*
* METHOD   : We first make the appropriate knot vectors, then we calulate
*            parameter values for the interpolation. The evaluator evalp
*            is used to calculate the interpolation points at these parameter
*            values. At last these points are interpolated.
*
*
* REFERENCES :
*
*-
* CALLS      : s2533(),evalp(),s1891().
*
* WRITTEN BY   :  Ulf J Krystad, SINTEF, Oslo, Norway.     Date: 1995-1
* REWRITTEN BY :  Johannes Kaasa, SINTEF, Oslo, Norway.    Date: 1995-8
*
*********************************************************************
*/
{
   int ki, kj;            /* Indices.                                   */
   int newin1;            /* Number of coefficents in first direction.  */
   int newin2;            /* Number of coefficents in second direction. */
   int leftknot1;         /* Pointer into knot array.                   */
   int leftknot2;         /* Pointer into knot array.                   */
   int open;              /* Open flag.                                 */
   int local_in;          /* Local number of coefficients.              */
   int nlr = 0;           /* Parameter to s1891.                        */
   int nrc = 0;           /* Parameter to s1891.                        */
   int *u_eder = SISL_NULL;    /* Parametrization of derivatives.            */
   int *v_eder = SISL_NULL;    /* Parametrization of derivatives.            */
   double par[2];         /* Surface parameters.                        */
   double *newet1 = SISL_NULL; /* Knot vector in first direction.            */
   double *newet2 = SISL_NULL; /* Knot vector in second direction.           */
   double *coef1 = SISL_NULL;  /* Surface coefficients.                      */
   double *coef2 = SISL_NULL;  /* Surface coefficients.                      */
   double *coef3 = SISL_NULL;  /* Surface coefficients.                      */
   double *u_par = SISL_NULL;  /* Schoenberg parameters in first direction.  */
   double *v_par = SISL_NULL;  /* Schoenberg parameters in second direction. */


   /* Check input */
   
   if (surf == SISL_NULL || u_multinc < 0 || v_multinc < 0 ||
       newik1 < (u_multinc + 2) || newik2 < (v_multinc + 2)) goto err150;

   /* Generate the knot array (spline space) in first direction. */

   s2533 (surf->et1, surf->ik1, surf->in1, u_multinc, newik1, &newin1, 
	  &newet1, stat);
   if (*stat < 0) goto error;

   /* Generate the knot array (spline space) in second direction. */

   s2533 (surf->et2, surf->ik2, surf->in2, v_multinc, newik2, &newin2, 
	  &newet2, stat);
   if (*stat < 0) goto error;
   
   /* Allocate utility arrays. */
   
   if ((coef1 = newarray(newin1*newin2*eval_dim, DOUBLE)) == SISL_NULL) 
      goto err101;
   if ((u_par = newarray(newin1, DOUBLE)) == SISL_NULL) goto err101;
   if ((v_par = newarray(newin2, DOUBLE)) == SISL_NULL) goto err101;
   if ((u_eder = newarray(newin1, INT)) == SISL_NULL) goto err101;
   if ((v_eder = newarray(newin2, INT)) == SISL_NULL) goto err101;
   
   /* Evaluate the property in the Schoenberg points. */
   
   for (ki = 0; ki < newin1; ki++)
   {
      u_par[ki] = 0.;
      for (kj = 1; kj < newik1; kj++)
	 u_par[ki] += newet1[ki + kj];
      u_par[ki] /= (newik1 - 1);
      
      u_eder[ki] = 0;
   }
   
   for (ki = 0; ki < newin2; ki++)
   {
      v_par[ki] = 0.;
      for (kj = 1; kj < newik2; kj++)
	 v_par[ki] += newet2[ki + kj];
      v_par[ki] /= (newik2 - 1);
      
      v_eder[ki] = 0;
      
      par[1] = v_par[ki];
      for (kj = 0; kj < newin1; kj++)
      {
	 par[0] = u_par[kj];
	 evalp(surf, 0, 1, 1, par, &leftknot1, &leftknot2,
	       &coef1[(ki*newin1 + kj)*eval_dim], stat);
	 if (*stat < 0 || *stat == 2) goto error;
      }
   }
   
   /* Interpolate curves in 1. parameter direction.  */
  
   open = SISL_CRV_OPEN;
   s1891(u_par, coef1, eval_dim, newin1, newin2, u_eder, open, newet1,
	 &coef2, &local_in, newik1, nlr, nrc, stat);
   if (*stat < 0) goto error;
  
   /* Interpolation in 2. parameter direction.                */
  
   s1891(v_par, coef2, newin1*eval_dim, newin2, 1, v_eder, open, newet2,
	 &coef3, &local_in, newik2, nlr, nrc, stat);
   if (*stat < 0) goto error;

   /* Create surface.  */
		 
   if ((*rsurf = newSurf(newin1, newin2, newik1, newik2, newet1, newet2,
			 coef3, 1, eval_dim, 2)) == SISL_NULL) goto err101;

  
   goto out;
  
  
  
   /* ---------------------- ERROR EXITS ------------------------------- */

   /* Error in space allocation */
   
 err101: 
   *stat = -101;
   s6err("s2534", *stat, 0);
   goto out;

   /* Error in input. */
   
 err150:
   *stat = -150;
   s6err("s2534", *stat, 0);
   goto out;

   /* Error in lower level routine. */
   
 error:
   s6err("s2534", *stat, 0);
   goto out;

   /* ---------------------- NORMAL EXIT ------------------------------- */

 out:
   if (coef1 != SISL_NULL)  freearray(coef1);
   if (coef2 != SISL_NULL)  freearray(coef2);
   if (u_par != SISL_NULL)  freearray(u_par);
   if (v_par != SISL_NULL)  freearray(v_par);
   if (u_eder != SISL_NULL) freearray(u_eder);
   if (v_eder != SISL_NULL) freearray(v_eder);
   
   return;
}
