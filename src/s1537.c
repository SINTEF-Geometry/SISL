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

#define S1537

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1537(double points[],int im1,int im2,int idim,double par1[],
	   double par2[],int con1,int con2,int con3,int con4, 
	   int order1, int order2,int iopen1,int iopen2,
	   SISLSurf **rsurf,int *jstat)
#else
void s1537(points,im1,im2,idim,par1, par2,con1,con2,con3,
	   con4,order1, order2,iopen1,iopen2,rsurf,jstat)
     double points[];
     int im1;
     int im2;
     int idim;
     double par1[];
     double par2[];
     int con1;
     int con2;
     int con3;
     int con4;
     int order1;
     int order2;
     int iopen1;
     int iopen2;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute a B-spline tensor surface interpolating a set
*          of points.
*
* INPUT:
*          points - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
*          par1   - Parametrization in first parameter direction.
*                   For closed curves, one additional parameter value
*                   must be spesified. The last entry contains
*                   the parametrization of the repeted start point.
*                   (if the endpoint is equal to the startpoint of
*                   the interpolation the lenght of the array could
*                   be equal to im1 also in the closed case).
*
*          par2   - Parametrization in second parameter direction.
*                   For closed curves, one additional parameter value
*                   must be spesified. The last entry contains
*                   the parametrization of the repeted start point.
*                   (if the endpoint is equal to the startpoint of
*                   the interpolation the lenght of the array could
*                   be equal to im2 also in the closed case).
*
*
*                          ^ Second par. direction 
*                          |     
*                          |    (2.)
*                          |-----------|
*                          |           |
*                     (3.) |           | (4.) 
*                          |           |
*                          |           |
*                          |-----------|-> First par. direction
*                               (1.)
*
*          con1      - Additional condition along edge 1:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con2      - Additional condition along edge 2:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con3      - Additional condition along edge 3:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con4      - Additional condition along edge 4:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          order1    - Order of surface in first parameter direction.
*
*          order2    - Order of surface in second parameter direction.
*
*          iopen1    - Open/close parameter in first parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed, non-periodic surface.
*                      = -1 : periodic surface.
*
*          iopen2    - Open/close parameter in second parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed, non-periodic surface.
*                      = -1 : periodic surface.
*
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     The interpolation is accomplished by using a one dimensional
*     routine for spline interpolation called several times. 
*     First, the datapoints
*     are considered to be idim*im1 dimentional and so on...
*
*
* REFERENCES :
*
* CALLS      : s1357
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, May 1993.
* REWISED BY : Vibeke Skytt, SINTEF, 0394. Introduced iopen1, iopen2.
*
*********************************************************************
*/                                                               
{
  int i;              /* Loop control parameter                      */
  int maxim;          /* Max (im1, im2)                              */
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  int newin1, newin2; /* Number of vertices along par. dir. 1 & 2    */
  int numpt;              /* Needed in call to s1357                 */
  double start=0;         /* Needed in call to s1357                 */
  double end;
  int *typept=SISL_NULL;       /* Array needed for call to s1357          */
  double *pointpar=SISL_NULL;  /* Array needed for call to s1357          */
  double *newcoeff=SISL_NULL;  /* Array needed for call to s1357          */  
  SISLCurve *curve1=SISL_NULL, *curve2=SISL_NULL;
  
  /* Allocate necessary array for call to s1357 */

  maxim = MAX( im1, im2 );
  if((typept = newarray(maxim, INT))==SISL_NULL) goto err101;
  for(i=0; i<maxim; i++)
    typept[i] = 1;

  /* Interpolation in second direction */

  s1357(points, im2, idim*im1, typept, par2, con1, con2, iopen2, order2, 
	start, &end, &curve1, &pointpar, &numpt, &kstat);
  if(kstat < 0) goto error;
  if(pointpar != SISL_NULL) 
    {
      freearray(pointpar);
      pointpar = SISL_NULL;
    }

  newin2 = curve1->in;

  /* Transpose result, store new coefficients in 
   * array newcoeff */

  if( (newcoeff = newarray(idim * im1 * newin2, DOUBLE)) == SISL_NULL )
    goto err101;
  s6chpar(curve1->ecoef, im1, newin2, idim, newcoeff);

  /* Interpolation in first parameter direction */

  s1357(newcoeff, im1, idim*newin2, typept, par1, con3, con4, iopen1, order1, 
	start, &end, &curve2, &pointpar, &numpt, &kstat);
  if(kstat < 0) goto error;
  if(pointpar != SISL_NULL) 
    {
      freearray(pointpar);
      pointpar = SISL_NULL;
    }

  newin1 = curve2->in;

  /* Transpose back coefficients */

  if( (newcoeff=increasearray(newcoeff, idim*newin1*newin2, DOUBLE)) 
     == SISL_NULL )  goto err101;
  s6chpar(curve2->ecoef, newin2, newin1, idim, newcoeff);

  /* Create instance of surface */

  if (((*rsurf) = newSurf(newin1, newin2, order1, order2, curve2->et,
		     curve1->et, newcoeff, 1, idim, 1)) == SISL_NULL)
     goto err101;
  
  /* Set periodicity flag.  */
  
  (*rsurf)->cuopen_1 = curve2->cuopen;
  (*rsurf)->cuopen_2 = curve1->cuopen;

  /* Success */
  
  *jstat = 0;
  goto out;  
  
  /* Allocation error. */

  err101: 
    *jstat = -101;
    s6err("s1537",*jstat,kpos);
    goto out;
      
  /* Error in lower level routine. */

  error:  *jstat =kstat;
    s6err("s1537",*jstat,kpos);
    goto out;
  
  out:
    /* Free arrays */
  
    if (newcoeff != SISL_NULL) freearray(newcoeff);
    if (typept != SISL_NULL) freearray(typept);

    /* Free local SISL-curve objects */
   
    if (curve1 != SISL_NULL) freeCurve(curve1);
    if (curve2 != SISL_NULL) freeCurve(curve2);
  
    return;
}
