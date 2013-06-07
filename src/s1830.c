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
 * $Id: s1830.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1830

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1830(SISLSurf *psurf,SISLCurve *pcurve,int *jstat)
#else
void s1830(psurf,pcurve,jstat)
     SISLSurf  *psurf;
     SISLCurve *pcurve;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform improved box-test in curve/surface intersection.
*              The test is based on the main tangent of the curve and
*              the main normal of the surface.
*
*
*
* INPUT      : psurf  - Pointer to surface involved in intersection problem.
*              pcurve - Pointer to curve involved in intersection problem.
*
*
*
* OUTPUT     : jstat  - status messages  
*                              = 2      : Only edge-intersections possible.
*                              = 1      : Possibility of intersection.
*                              = 0      : No possibility of intersection.
*                              < 0      : error
*
*
* METHOD     : Let the main tangent of the curve be the vector from the
*              first to the last vertex, and let the main normal of the
*              surface be the cross-product of the main diagonals. Perform
*              rotated box-tests with respect to this two vectors.
*
*
* REFERENCES :
*
*-
* CALLS      : s1834  - Perform rotated box-test.
*              s6diff - Difference vector between two vectors.
*              s6crss - Cross-product of two vectors.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  int kstat = 0;         /* Local status variable.           */
  int kpos = 0;          /* Position of error.               */
  int kdim;              /* Dimension of space.              */
  int knc;               /* Number of vertices of curve.     */
  int kn1,kn2;           /* Number of vertices of surface.   */
  double *scurve;        /* Vertices of curve.               */
  double *ssurf;         /* Vertices of surface.             */
  double *stan = SISL_NULL;   /* Main tangent of curve.           */
  double *sdiag1 = SISL_NULL; /* First main diagonal of surface.  */
  double *sdiag2 = SISL_NULL; /* Second main diagonal of surface. */
  double *snorm = SISL_NULL;  /* Main normal of surface.          */
  
  /* Test input.  */
  
  kdim = psurf -> idim;
  if (kdim != 3) goto err104;
  if (kdim != pcurve -> idim) goto err106;
  
  /* Allocate space for local arrays.  */
  
  if ((stan = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((sdiag1 = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((sdiag2 = newarray(kdim,double)) == SISL_NULL) goto err101;
  if ((snorm = newarray(kdim,double)) == SISL_NULL) goto err101;
  
  /* Describe curve with local parameters.  */
  
  knc = pcurve->in;
  scurve = pcurve->ecoef;
  
  /* Describe surface with local parameters.  */
  
  kn1 = psurf->in1;
  kn2 = psurf->in2;
  ssurf = psurf->ecoef;
  
  /* Fetch main tangent of curve.  */
  
  s6diff(scurve+(knc-1)*kdim,scurve,kdim,stan);
  
  /* Fetch main diagonals of surface.  */
  
  s6diff(ssurf+(kn1*kn2-1)*kdim,ssurf,kdim,sdiag1);
  
  s6diff(ssurf+kn1*(kn2-1)*kdim,ssurf+(kn1-1)*kdim,kdim,sdiag2);
  
  /* Compute main normal of surface.  */
  
  s6crss(sdiag1,sdiag2,snorm);
  
  /* Perform rotated box-test.  */
  
  s1834(ssurf,kn1*kn2,scurve,knc,kdim,stan,snorm,&kstat);
  if (kstat < 0) goto error;
  
  if (kstat == 1)
    {
      kstat = 0;
      
      s1834(ssurf,kn1*kn2,scurve,knc,kdim,snorm,stan,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* Improved box-test performed.  */
  
  *jstat = kstat;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1830",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 3.  */
  
 err104: *jstat = -104;
  s6err("s1830",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("s1830",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1830",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free space occupied by local arrays.  */
  
  if (stan != SISL_NULL) freearray(stan);
  if (sdiag1 != SISL_NULL) freearray(sdiag1);
  if (sdiag2 != SISL_NULL) freearray(sdiag2);
  if (snorm != SISL_NULL) freearray(snorm);
  
  return;
}
