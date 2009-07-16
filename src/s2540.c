/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1995 by                                                     */
/*     SINTEF, Oslo, Norway.                                                 */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s2540.c,v 1.4 1999-01-06 12:23:42 jka Exp $
 *
 */

#define S2540

#include "sislP.h"



#if defined(SISLNEEDPROTOTYPES)
void
s2540(SISLSurf *surf, int curvature_type, int export_par_val, int pick_subpart,
      double boundary[], int n_u, int n_v,
      double **garr, int *stat)
#else
void
s2540(surf, curvature_type, export_par_val, pick_subpart, boundary, n_u, n_v,
      garr, stat)
     SISLSurf *surf;
     int curvature_type;
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
* PURPOSE : To compute a set of curvature values on an uniform grid
*           in a selected subset of the parameter domain for a
*           NURBS surface.
*
*
* INPUT   : surf  	   - The surface to evaluate.
*	    curvature_type - The type of curvature:
*                            0 - Gaussian.
*                            1 - Mean.
*                            2 - Absolute.
*                            3 - Total.
*                            4 - second order Mehlum (curvature).
*                            5 - third order Mehlum (variation of curvature).
*	    export_par_val - Flag telling if the parameter values for each grid
*                            point is to be exported:
*                            0 - False, do not export parameter values,
*                            1 - True, do export parameter values.
*           pick_subpart   - Flag teliing if the grid is to be calculated on a
*                            subpart of the surface:
*                            0 - False, calculate grid on the complete surface,
*                            1 - True, calculate grid on a part of the surface.
*           boundary       - A rectangular subset of the parameter domain.
*              	             [0] - Min value 1. parameter direction.
*                            [1] - Min value 2. parameter direction.
*                            [2] - Max value 1. parameter direction.
*                            [3] - Max value 2. parameter direction.
*                            ONLY USED WHEN pick_subpart = 1. If pick_subpart
*                            = 0, the parameter area of surf is given out here.
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
**
*           stat           - Status message.
*                               > 0      : Warning.
*                               = 0      : Ok.
*                               < 0      : Error.
*
*
* METHOD  :
*
*
* CALLS   :s2541(), ( s2500(), s2502(), s2504(), s2506(), s2508(), s2510()).
*
* WRITTEN BY :  Ulf J Krystad,  SINTEF, Oslo, Norway, Jan. 1995.
* REVISED BY :  Johannes Kaasa, SINTEF, Oslo, Norway, Aug. 1995.
*               (Mehlum curvature added).
*
*********************************************************************
*/
{
  int dim=1;              /* Dimension of evaluator output.	  */

  if (curvature_type == 0)
     s2541(surf, s2500, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else if (curvature_type == 1)
     s2541(surf, s2502, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else if (curvature_type == 2)
     s2541(surf, s2504, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else if (curvature_type == 3)
     s2541(surf, s2506, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else if (curvature_type == 4)
     s2541(surf, s2508, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else if (curvature_type == 5)
     s2541(surf, s2510, dim, export_par_val, pick_subpart, boundary, n_u, n_v,
	   garr, stat);
  else
     goto err151;
  if (*stat < 0) goto error;


  *stat = 0;
  goto out;

  /* ___________________________________________________________________ */
  /*                         ERROR EXITS                                 */
  /* ___________________________________________________________________ */

  /* Error in input, uknown type specified */
err151:
  *stat = -151;
  s6err("s2540", *stat, 0);
  goto out;

  /* Error in lower level routine. */
error:
  s6err("s2540", *stat, 0);
  goto out;

  /* ___________________________________________________________________ */
  /*                         THE ONE AND ONLY EXIT                       */
  /* ___________________________________________________________________ */

out:
  return;

}
