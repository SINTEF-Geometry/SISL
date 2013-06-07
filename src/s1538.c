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

#define S1538

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s1538(int inbcrv,SISLCurve *vpcurv[],int nctyp[],double astpar,
	   int iopen,int iord2,int iflag,
	   SISLSurf **rsurf,double **gpar,int *jstat)
#else
void s1538(inbcrv,vpcurv,nctyp,astpar,iopen,iord2,
           iflag,rsurf,gpar,jstat)
     int    	inbcrv;
     SISLCurve  *vpcurv[];
     int   	nctyp[];
     double	astpar;
     int    	iopen;
     int    	iord2;
     int    	iflag;
     SISLSurf   **rsurf;
     double 	**gpar;
     int    	*jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To create a spline lofted surface
*              from a set of input-curves.
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcurv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*              nctyp  - Array (length inbcrv) containing the types
*                       of curves in the curve-set.
*                        1 - Ordinary curve.
*                        2 - Knuckle curve. Treated as ordinary curve.
*                        3 - Tangent to next curve.
*                        4 - Tangent to prior curve.
*                       (5 - Double derivative to prior curve.)
*                       (6 - Double derivative to next curve.)
*                       13 - SISLCurve giving start of tangent to next curve.
*                       14 - SISLCurve giving end of tengent to prior curve.
*              astpar - Start-parameter for spline lofting direction.
*              iopen  - Flag telling if the resulting surface should
*                       be closed or open.
*                       -1 - The surface should be closed and periodic.
*                        0 - The surface should be closed.
*                        1 - The surface should be open.
*              iord2  - Maximal order of the B-spline basis in the
*                       lofting direction.
*              iflag  - Flag telling if the size of the tangents in the
*                       derivative curves should be adjusted or not.
*                        0 - Do not adjust tangent-sizes.
*                        1 - Adjust tangent-sizes.
*
* OUTPUT     : rsurf  - Pointer to the surface produced.
*              gpar   - The input-curves are constant parameter-lines
*                       in the parameter-plane of the produced surface.
*                       (i) - contains the (constant) value of this
*                             parameter of input-curve no. i.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : A common basis for all the B-spline curves are found.
*              The curves are represented using this basis.
*              The resulting curves are given to an interpolation
*              routine that calculates the B-spline vertices of the
*              resulting spline lofted surface.
*              Throughout these routines, first parameterdirection
*              will be the interpolating direction, second parameter-
*              direction will be along the input curves.
*-
* CALLS      : s1931,s1917,s1918,s1358,s6err.
*
* WRITTEN BY : A. M. Ytrehus   SI  Oslo,Norway. Sep. 1988
* Revised by : Tor Dokken, SI, Oslo, Norway, 26-feb-1989
* Revised by : Trond Vidar Stensby, SI, 91-08
* REVISED BY: Vibeke Skytt, 03.94. This routine corresponds to s1333,
*                                  but differ in the use of the parameter
*                                  iopen.
* Revised by : Paal Fugelli, 17/08-1994.  Fixed memory leak from 'gpar'
*              allocated in s1357().
*
*********************************************************************
*/
{
  int kind,kcopy,kdim;
  int kn1,kord1,knbcrv;
  int kcnsta,kcnend;         /* Interpolation condition at start or end */
  int ki,kj,kl,km;
  int kleng;                 /* Number of doubles describing a curve  */
  int ktype;                 /* Kind of interpolation condition.      */
  int kopen;                 /* Open/closed parameter in curve direction. */
  SISLCurve *qc;             /* Pointer to curve representing surface */
  int *lder = SISL_NULL;	     /* Derivative indicators from s1915. */
  double *spar=SISL_NULL; 	     /* Param. values of point conditions. */
  double *spar2=SISL_NULL; 	     /* Parameter values from s1915. */
  double *sknot1=SISL_NULL;       /* Knot vector.                 */
  double *scoef2=SISL_NULL;       /* Pointer to vertices expressed in same basis  */
  double tstpar;             /* Parameter value of last curve                */
  int kstat = 0;             /* Status variable. */
  int kpos = 0;              /* Position of error. */
  int knbpar;                /* Number of parameter values produced          */
  int kdimcrv;               /* kdim multiplied with number of vertices kn1  */
  int kcont;                 /* Continuity at end of curves */


  /* Initiate variables. */

  kdim = vpcurv[0]->idim;

  if (inbcrv < 2) goto err179;

  /* Put the curves into common basis. */

  s1931 (inbcrv, vpcurv, &sknot1, &scoef2, &kn1, &kord1, &kstat);
  if (kstat < 0)
    goto error;

  /* Create the parameter-values for the knot-vector
     (in lofting direction) for a lofted surface, allocate array for
     parameter values.    */

  s1917 (inbcrv, scoef2, kn1, kdim, nctyp, astpar, iopen,
	 &spar2, &lder, &knbcrv, &kstat);

  if (kstat < 0)
    goto error;

  /* Convert condition 13 and 14 to 3 and 4 */

  kleng = kn1*kdim;
  for (ki=0 ; ki<knbcrv ; ki++)
    {
       ktype = nctyp[ki];

      if (ktype == 13 && ki+1<knbcrv)
        {
	  /*
	   * Start of tangent to next curve,
	   * make difference of next curve and this curve
	   */

	  for (kj=ki*kleng,kl=kj+kleng,km=0; km <kleng ; kj++,kl++,km++)
	      scoef2[kj] = scoef2[kl] - scoef2[kj];
	  nctyp[ki] = 3;
        }
      else if (ktype == 14 && ki>0)
        {
	  /* End of tangent to prior curve,
	   * make difference of this curve
	   * and prior curve
	   */

	  for (kj=ki*kleng,kl=kj-kleng,km=0; km <kleng ; kj++,kl++,km++)
	      scoef2[kj] = scoef2[kj] - scoef2[kl];
	  nctyp[ki] = 4;
        }
    }

  spar = newarray(knbcrv+1,DOUBLE);
  if (spar==SISL_NULL) goto err101;

  /*  Only copy parameter values of point conditions */

  for (ki=0,kl=0; ki<knbcrv ; ki++)
    {
      if (nctyp[ki] == 1 || nctyp[ki] == 2)
        {
	  spar[kl] = spar2[ki];
	  kl++;
        }
    }

  /* Add one extra parameter value if closed curve */

  if (iopen != SISL_CRV_OPEN) spar[kl] = spar2[knbcrv];

  /* Adjust tangent-lengths if wanted. */

  if (iflag)
    {
      s1918 (knbcrv, sknot1, scoef2, kn1, kord1, kdim, spar2, lder, &kstat);
      if (kstat < 0) goto error;
    }

  /* Interpolate with point interpolation method */

  kcnsta = 0;
  kcnend = 0;
  kdimcrv = kdim*kn1;

  s1357(scoef2,knbcrv,kdimcrv,nctyp,spar,kcnsta,kcnend,iopen,iord2,astpar,
	&tstpar,&qc,gpar,&knbpar,&kstat);
  if (kstat<0) goto error;

  /* The knot vector in the lofting direction and the coefficients are
     now contained in the curve object pointed to by qc */

  /* Create the surface */

  kind = 1;
  kcopy = 1;
  *rsurf = newSurf(kn1,qc->in,kord1,qc->ik,sknot1,qc->et,qc->ecoef,
		   kind,kdim,kcopy);
  if (*rsurf == SISL_NULL) goto err101;

  /* Copy cuopen flag from curve */
  (*rsurf)->cuopen_2 = qc->cuopen;

  /* Release the curve object */

  freeCurve(qc);

  /* Output parametervalues according to the input curves, but must
     remember to free the space allocated in call to s1357() first.  */

  if ( (*gpar) != SISL_NULL ) freearray(*gpar);  /* PFU 17/08-94. */
  *gpar = spar;

  /* Decide if the surface should have a cyclic behaviour in first
     parameter direction i.e. the direction of the curves */

  s1333_count(inbcrv,vpcurv,&kcont,&kstat);
  if (kstat<0) goto error;

  if (kcont>=0)
      {
        s1333_cyclic(*rsurf,kcont,&kstat);
	if (kstat<0) goto error;

	/* Set periodic flag */
	(*rsurf)->cuopen_1 = SISL_SURF_PERIODIC;
      }
      else
      {
         /* Test if the surface should be closed and non-periodic.  */

         for (kopen=-2, ki=0; ki<inbcrv; ki++)
           kopen = MAX(kopen,vpcurv[ki]->cuopen);
         if (kopen == SISL_CRV_CLOSED) (*rsurf)->cuopen_1 = SISL_SURF_CLOSED;
      }

  /* Task done */

  *jstat = 0;
  goto out;

  /* Error in allocation. */

 err101:
  *jstat = -101;
  s6err("s1538",*jstat,kpos);
  goto out;


  /* Error in interpolation conditions. No. of curves < 2. */

 err179:
  *jstat = -179;
  s6err("s1538",*jstat,kpos);
  goto out;


  /* Error in lower level routine.  */

  error :
    *jstat = kstat;
  s6err("s1538",*jstat,kpos);
  goto out;
 out:

  /* Free allocated scratch  */

  if (sknot1 != SISL_NULL) freearray(sknot1);
  if (scoef2 != SISL_NULL) freearray(scoef2);
  if (spar2 != SISL_NULL) freearray(spar2);
  if (lder != SISL_NULL) freearray(lder);

  return;
}
