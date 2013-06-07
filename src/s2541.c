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
 * $Id: s2541.c,v 1.4 2005-02-28 09:04:49 afr Exp $
 *
 */

#define S2541

#include "sislP.h"




#if defined(SISLNEEDPROTOTYPES)
void
   s2541(SISLSurf *surf,
	 void evalp(SISLSurf *surf, int ider,int iside1,int iside2,
		    double parvalue[], int *leftknot1, int *leftknot2,
		    double *result, int *istat),
	 int dim,
	 int export_par_val,
	 int pick_subpart,
	 double boundary[],
	 int n_u,
	 int n_v,
	 double **garr,
	 int *stat)
#else
void
   s2541(surf, evalp, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	 garr, stat)
      SISLSurf *surf;
      void evalp();
      int dim;
      int export_par_val;
      int pick_subpart;
      double boundary[];
      int n_u;
      int n_v;
      double **garr;
      int *stat;
#endif
/*
***************************************************************************
*
***************************************************************************
* PURPOSE : To compute a set of values on an uniform grid
*           in a selected subset of the parameter domain for a
*           NURBS surface. The values given is based on the evaluator evalp().
*
*
* INPUT   : surf           - The surface to evaluate.
*           evalp   	   - A pointer to the function/evaluator to be used
*                            when computing the wanted values.
*           dim    	   - The spatial dimension required for the output
*                            values, must be consistent with the result from
*                            the evaluator.
*	    export_par_val - Flag telling if the parameter values for each grid
*                            point is to be exported:
*                            0 - False, do not export parameter values,
*                            1 - True, do export parameter values.
*           pick_subpart   - Flag teliing if the grid is to be calculated on a
*                            subpart of the surface:
*                            0 - False, calculate grid on the complete surface,
*                            1 - True, calculate grid on a part of the surface.
*           boundary       - A rectangular subset of the parameter domain.
*              		     [0] - Min value 1. parameter direction.
*                            [1] - Min value 2. parameter direction.
*                            [2] - Max value 1. parameter direction.
*                            [3] - Max value 2. parameter direction.
*                            ONLY USED WHEN pick_subpart = 1.
*           n_u            - Number of segments in 1. parameter direction.
*           n_v            - Number of segments in 2. parameter direction.
*
*
*
* OUTPUT  : garr  	   - Array containing the computed values on the grid.
*		             The allocation is done internally and the dimension
*			     is  3*(n_u+1)*(n_v+1) if export_par_val is true,
*		             and (n_u+1)*(n_v+1) if export_par_val is false.
*                            Each gridpoint consist of a triple
*                            (Ui,Vj,curvature(Ui,Vj)) or only curvature(Ui,Vj).
*			     The sequence is running first in the
*                            1. parameter direction.
*
*           stat           - Status message.
*                               > 0      : Warning.
*                               = 0      : Ok.
*                               < 0      : Error.
*
*
* METHOD  :
*
*
* CALLS   :s1001(), ( s2500(), s2502(), s2504(), s2506(), s2508()).
*
* WRITTEN BY :  Ulf J Krystad, SINTEF, Oslo, Norway, Jan. 1995.
* REVISED BY :  Johannes Kaasa, SINTEF, Oslo, Norway, Aug. 1995.
*
*********************************************************************
*/
{
  int dimpnt;			/* Dimesion of each gridpoint        	*/
  int kder = 0;                 /* Derivativ indicator          	*/
  int ki, kj;               	/* Loop control variable        	*/
  int klfs = 0;                 /* Pointer into knot vector           	*/
  int klft = 0;            	/* Pointer into knot vector   		*/
  int kside1=0,kside2=0;    	/* Left,right evaluations             	*/
  int incr;    			/* Increment, space for par values.    	*/
  double duv[2];		/* Increment in 1.+2. par dir.  	*/
  double UV[2];			/* Current grid value in par space.	*/
  double *sarr = SISL_NULL;    	/* Local pointer eq (*garr)		*/
  double *sp = SISL_NULL;      	/* Local pointer into (*garr)		*/
  SISLSurf *temp = SISL_NULL;  	/* Temp surface. 			*/
  /* __________________________________________________________________ */


  /* Initiate output variables . */

  *garr  = SISL_NULL;
  *stat = 0;


  /* Check input. */
  if ( !surf ) 		goto err150;
  if (dim < 1) 	        goto err102;
  if ( n_u < 1) 	goto err172;
  if ( n_v < 1) 	goto err172;

  /* Pick the surface defined over the subset wanted. */

  if (pick_subpart == 1)
  {
     s1001 (surf, boundary[0], boundary[1], boundary[2], boundary[3],
   	 &temp, stat);
     if (*stat < 0) goto error;
  }
  else
  {
     temp = surf;
     boundary[0] = temp->et1[temp->ik1 - 1];
     boundary[1] = temp->et2[temp->ik2 - 1];
     boundary[2] = temp->et1[temp->in1];
     boundary[3] = temp->et2[temp->in2];
  }

  /* Allocate space needed */
  incr = (export_par_val? 2 : 0);
  dimpnt = incr+dim;
  if ((sarr = newarray(dimpnt*(n_u+1)*(n_v+1), double)) == SISL_NULL) goto err101;

  /* The evaluation loop, note that to ensure that we get boundary values
     correct, we stop the main loops one step too early and then jump
     to the maximum boundary. */

  duv[0] = (boundary[2] - boundary[0])/n_u;
  duv[1] = (boundary[3] - boundary[1])/n_v;

  for (kj=0,UV[1]=boundary[1],sp=sarr;kj<n_v;kj++,UV[1] += duv[1])
  {
     for (ki=0,UV[0]=boundary[0];ki<n_u;ki++,UV[0] += duv[0],sp += dimpnt)
     {
	if (export_par_val)
	{
	   sp[0] = UV[0];
	   sp[1] = UV[1];
	}
	evalp(temp,kder=0,kside1=0,kside2=0,UV,&klfs,&klft,
	      sp+incr, stat);
	if (*stat < 0) goto error;
     }

     /* Last column. */

     UV[0] = boundary[2];
     if (export_par_val)
     {
	sp[0] = UV[0];
	sp[1] = UV[1];
     }
     evalp(temp,kder=0,kside1=0,kside2=0,UV,&klfs,&klft,
	   sp+incr, stat);
     if (*stat < 0) goto error;
     sp += dimpnt;
  }

  /* Last row. */

  for (ki=0,UV[0]=boundary[0];ki<n_u;ki++,UV[0] += duv[0],sp += dimpnt)
  {
     UV[1] = boundary[3];
     if (export_par_val)
     {
	sp[0] = UV[0];
	sp[1] = UV[1];
     }
     evalp(temp,kder=0,kside1=0,kside2=0,UV,&klfs,&klft,
	   sp+incr, stat);
     if (*stat < 0) goto error;
  }

  /* Last column in last row. */

  UV[0] = boundary[2];
  UV[1] = boundary[3];
  if (export_par_val)
  {
     sp[0] = UV[0];
     sp[1] = UV[1];
  }
  evalp(temp,kder=0,kside1=0,kside2=0,sp,&klfs,&klft,
	sp+incr, stat);
  if (*stat < 0) goto error;


  /* OK, we are thru, hand over the array: */
  *garr  = sarr;
  sarr   = SISL_NULL;
  *stat = 0;

  goto out;



  /* ___________________________________________________________________ */
  /*                         ERROR EXITS                                 */
  /* ___________________________________________________________________ */


  /* Error in space allocation */
err101:
  *stat = -101;
  s6err("s2541", *stat, 0);
  goto out;

   /* Error in input, dim < 1 */
err102:
   *stat = -102;
   s6err("s2541", *stat, 0);
   goto out;

   /* Error. Input (surface) pointer is SISL_NULL. */
err150:
  *stat = -150;
  s6err("s2541", *stat, 0);
  goto out;

  /* Error. Too few segments in input. */
err172:
  *stat = -172;
  s6err("s2541", *stat, 0);
  goto out;

  /* Error in lower level routine. */
error:
  s6err("s2541", *stat, 0);
  goto out;

  /* ___________________________________________________________________ */
  /*                         THE ONE AND ONLY EXIT                       */
  /* ___________________________________________________________________ */

out:
  /* Free local heap space. */

  if (pick_subpart == 1 && temp) freeSurf(temp);
  if ( sarr )  		         freearray(sarr);
  return;

}
