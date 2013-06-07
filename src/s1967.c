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

#define S1967

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1967(double ep[],double etang1[],double etang2[],double eder11[],
	   int im1,int im2,int idim,int ipar,double epar1[],double epar2[],
	   double eeps[],int nend[],int iopen1,int iopen2,double edgeps[],
	   int iopt,int itmax,SISLSurf **rs,double emxerr[],
	   int *jstat)
#else
void s1967(ep,etang1,etang2,eder11,im1,im2,idim,ipar,epar1,epar2,
	   eeps,nend,iopen1,iopen2,edgeps,iopt,itmax,rs,emxerr,jstat)
     double ep[];
     double etang1[];
     double etang2[];
     double eder11[];
     int    im1;
     int    im2;
     int    idim;
     int    ipar;
     double epar1[];
     double epar2[];
     double eeps[];
     int    nend[];
     int    iopen1;
     int    iopen2;
     double edgeps[];
     int    iopt;
     int    itmax;
     SISLSurf   **rs;
     double emxerr[];
     int    *jstat;
#endif
/*
********************************************************************
*
* Purpose: To compute a bicubic hermite spline-approximation to the
*          position and derivative data given by ep,etang1,etang2
*          and eder11.
*
* Input : Ep     - Array (length idim*im1*im2) containing the points
*                  to be approximated.
*         Etang1 - Array (length idim*im1*im2) containing the derivatives
*                  (tangents) in the first parameter-direction at the
*                  data-points.
*         Etang2 - Array (length idim*im1*im2) containing the derivatives
*                  (tangents) in the second parameter-direction at the
*                  data-points.
*         Eder11 - Array (length idim*im1*im2) containing the cross (twist)
*                  derivatives at the data-points.
*         Im1    - The no. of points in the first parameter-direction.
*         Im2    - The no. of points in the second parameter-direction.
*         Ipar   - Flag determining the parametrization of the data-points.
*                   = 1: Mean accumulated cord length parametrization.
*                   = 2: Uniform parametrization.
*                   = 3: Parametrization given by epar1 and epar2.
*         Epar1  - Array (length im1) containing a parametrization in the
*                  first parameter-direction. (Will only be used if ipar=3.)
*         Epar2  - Array (length im2) containing a parametrization in the
*                  surface lies.)
*         Eeps   - Array (length idim) containing the maximum deviation
*                  which is acceptable in each of the idim components of
*                  the surface (except possibly along the edges).
*         Nend   - Array (length 4) containing the no. of derivatives 
*                  to be kept fixed along each edge of the surface.
*                  The numbering of the edges is the same as for edeps below.
*                  All the derivatives of order < nend(i)-1 will be kept fixed
*                  along edge no. i. Hence nend(i)=0 indicates that nothing
*                  is to be kpet fixed along edge no. i.
*                  To be kept fixed here means to have error less than edgeps.
*                  In general, it is impossible to remove any knots and
*                  keep an edge completely fixed.
*         iopen1 - Open/closed parameter in first parameter direction.
*                      =  1 : Produce open surface.
*                      =  0 : Produce closed, non-periodic surface if possible.
*                      = -1 : Produce closed, periodic surface if possible.
*                  NB! The surface will be closed/periodic only if the first 
*                      and last column of data are (approximately) equal.
*         iopen2 - Open/closed parameter in second parameter direction.
*                      =  1 : Produce open surface.
*                      =  0 : Produce closed, non-periodic surface if possible.
*                      = -1 : Produce closed, periodic surface if possible.
*                  NB! The surface will be closed/periodic only if the first 
*                      and last row of data are (approximately) equal.
*         Edgeps - Array (length idim*4) containing the max. deviation
*                  which is acceptable along the edges of the surfaces.
*                  Edgeps(1,i):edgeps(idim,i) gives the tolerance along
*                  the edge corresponding to the i-th parameter having
*                  its min. or max value.
*                  i=1 : min value of first parameter.
*                  i=2 : max value of first parameter.
*                  i=3 : min value of second parameter.
*                  i=4 : max value of second parameter.
*                  Edgeps(kp,i) will only have any significance if nend(i)>0.
*         Iopt   - Flag indicationg the order in which tha data-reduction
*                  is to be performed.
*                   = 1 : Remove knots in parameter-direction 1 only.
*                   = 2 : Remove knots in parameter-direction 2 only.
*                   = 3 : Remove knots in parameter-direction 1 and
*                         and then in parameter-direction 2.
*                   = 4 : Remove knots in parameter-direction 2 and
*                         and then in parameter-direction 1.
*         Itmax  - Max. no. of iteration.
*
* Ouput:  Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning.
*         Rs     - Pointer to surface.
*         Emxerr - Array (length idim) (allocated outside this routine.)
*                  containing an upper bound for the error comitted in
*                  each component during the data reduction.
*
* Method: First the bicubic hermite spline-interpolant is computed
*         using the appropriate parametrization, and then knots
*         are removed from this approximation by a call to the
*         data-reduction routine for surfaces.
*
* Calls: s1530, s1965, s6err.
*
* Written by: C.R.Birkeland, Si, Oslo, Norway, May 1993.
* Changed and renamed by : Vibeke Skytt, SINTEF Oslo, 02.95. Introduced
*                                                            periodicity.
**********************************************************************
*/
{
  int stat=0, kpos=0;         /* Error message parameters        */
  double *par1 = SISL_NULL;        /* Used to store parametrizations  */
  double *par2 = SISL_NULL;
  SISLSurf *osurf = SISL_NULL;     /* Hermite interp. surface         */ 

  /* Check Input */

  if (im1 < 2 || im2 < 2 || idim < 1) 
    goto err103; 
  if (ipar < 1 || ipar > 3) ipar = 1;

  /* Generate parametrization */

  if (ipar != 3) 
    { 
      /* Generate parametrization */
      
      s1528(idim, im1, im2, ep, ipar, SISL_CRV_OPEN, SISL_CRV_OPEN,
	    &par1, &par2, &stat);
      if (stat<0) goto error;
    }
  else 
    {
      /* Parametrization is passed as parameter */

      par1 = epar1;
      par2 = epar2;
    } 

  /* Perform bicubic hermite spline interpolation */

  s1530(ep, etang1, etang2, eder11, par1, par2,
	im1, im2, idim, &osurf, &stat);
  if (stat<0) goto error;

  /* Perform final datareduction step */

  s1965(osurf, eeps, nend, iopen1, iopen2, edgeps, iopt, itmax,  
        rs, emxerr, &stat);
  if (stat<0) goto error;
  
  /* Success */

  *jstat = 0;
  goto out;
  
  
  /* Error in input */

  err103: 
    *jstat = -103;
    s6err("s1967",*jstat,kpos);
    goto out;
  
  /* Error in lower level routine. */

  error: 
    *jstat = stat;
    s6err("s1967",*jstat,kpos);
    goto out;
  
  /* Exit. */

  out:
    /* Free SISL-surface allocated in this routine */

    if(osurf != SISL_NULL) freeSurf(osurf);

    if (ipar != 3)
      {
	if(par1 != SISL_NULL) freearray(par1);
	if(par2 != SISL_NULL) freearray(par2);
      }
    return;
}
