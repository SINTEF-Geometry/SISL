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
 * $Id:
 *
 */


#define S1508

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1508(int inbcrv, SISLCurve **vpcurv, double par_arr[],
       SISLSurf **rsurf, int *jstat)
#else
void
s1508(inbcrv, vpcurv, par_arr,
      rsurf, jstat)
     int inbcrv;
     SISLCurve **vpcurv;
     double par_arr[];
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To create a rational B-spline lofted surface
*              from a set of rational B-spline input-curves.
*
* INPUT      : inbcrv - Number of curves in the curve-set.
*              vpcurv  - Array (length inbcrv) of pointers to the
*                       curves in the curve-set.
*              par_arr - The required parameterisation, must be
*                        strictly increasing, length inbcrv.
*
* OUTPUT     : rsurf  - Pointer to the surface produced.
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error.
*
* METHOD     : A common basis for all the B-spline curves are found.
*              The curves are represented using this basis.
*              The resulting curves are given to an interpolation
*              routine that calculates the B-spline vertices of the
*              resulting spline lofted surface.
*              The surface will be C1 and cubic in the lofting direction.
*              Throughout these routines, first parameter direction
*              will be the interpolating direction, second parameter-
*              direction will be along the input curves.
*
* REFERENCES :
*
*-
* CALLS      : s1931,s1917,s1918,s1358,s6err.
*
* WRITTEN BY : Michael Floater, SI, Oslo, October 1993.
*
*********************************************************************
*/
{
  int kind, kcopy, kdim;
  int kn1, kord1;
  int i,ki,kj,kk;         /* Loop variables. */
  SISLCurve *qc;	  /* Pointer to curve representing surface */
  double *sknot1 = SISL_NULL;  /* Knot vector. */
  double *scoef2 = SISL_NULL;  /* Pointer to vertices expressed in same basis  */
  int kstat = 0;	  /* Status variable. */
  int kpos = 0;		  /* Position of error. */
  double *epcoef;
  double *evcoef;
  double *epweight;
  double *evweight;
  double *epweight2;
  double *evweight2;
  double *evnew;
  double *hvcoef;

  SISLCurve **ratcurves=SISL_NULL;



  *jstat = 0;


  /* Check input. */

  if(inbcrv < 2) goto err179;

  /* Initiate variables. */

  kdim = vpcurv[0]->idim;

  /* Convert the rational curves to homogeneous form.
     The temporary homogeneous curves only point to
     existing arrays so no major copying is required. */

  ratcurves = newarray(inbcrv,SISLCurve*);
  for(i=0; i<inbcrv; i++)
  {
    /* Just set pointers to the homogeneous coordinates. */
    ratcurves[i] = newCurve(vpcurv[i]->in,vpcurv[i]->ik,vpcurv[i]->et,
                            vpcurv[i]->rcoef,1,kdim+1,0);
  }

  /* Put the curves into common basis. */

  s1931 (inbcrv, ratcurves, &sknot1, &scoef2, &kn1, &kord1, &kstat);
  if (kstat < 0)
    goto error;

  /* scoef2 now holds homogeneous coordinates.
     Calculate control points and weights. */

  /* Allocate array for points and weights */

  epcoef    = newarray(kdim*inbcrv*kn1,DOUBLE);
  if (epcoef == SISL_NULL) goto err101;

  evcoef    = newarray(kdim*inbcrv*kn1,DOUBLE);
  if (evcoef == SISL_NULL) goto err101;

  epweight    = newarray(inbcrv*kn1,DOUBLE);
  if (epweight == SISL_NULL) goto err101;

  evweight    = newarray(inbcrv*kn1,DOUBLE);
  if (evweight == SISL_NULL) goto err101;

  epweight2    = newarray(inbcrv*kn1,DOUBLE);
  if (epweight2 == SISL_NULL) goto err101;

  evweight2    = newarray(inbcrv*kn1,DOUBLE);
  if (evweight2 == SISL_NULL) goto err101;


  for(ki=0,kj=0,kk=0; ki<inbcrv*kn1; ki++,kj+=(kdim+1),kk+=kdim)
  {
      for(i=0; i<kdim; i++)
      {
          epcoef[kk+i] = scoef2[kj+i] / scoef2[kj+kdim];
      }

      epweight[ki] = scoef2[kj+kdim];
  }


  /* Estimate derivative for control points and weights. */

  s1516(epcoef,par_arr,inbcrv,kn1*kdim,&evcoef,&kstat);
  if (kstat < 0) goto error;

  s1516(epweight,par_arr,inbcrv,kn1,&evweight,&kstat);
  if (kstat < 0) goto error;

  /* Adjust weight derivatives so that weight function is positive. */

  /* First change parameter direction. */

  s6chpar(epweight,kn1,inbcrv,1,epweight2);
  s6chpar(evweight,kn1,inbcrv,1,evweight2);

  for(ki=0,kk=0; ki<kn1; ki++,kk+=inbcrv)
  {
      s1517(epweight2+kk,evweight2+kk,par_arr,inbcrv,0.0,&evnew,&kstat);
      if (kstat < 0) goto error;

      memcopy(evweight2+kk,evnew,inbcrv,double);
      freearray(evnew);

  }

  /* Change back parameter direction. */

  s6chpar(evweight2,inbcrv,kn1,1,evweight);

  /* Allocate array for homogeneous derivatives. */

  hvcoef    = newarray(inbcrv*kn1*(kdim+1),DOUBLE);
  if (hvcoef == SISL_NULL) goto err101;

  /* Use product rule to estimate derivatives for
     homogeneous coordinates. */

  for(ki=0,kj=0,kk=0; ki<kn1*inbcrv; ki++,kj+=kdim,kk+=(kdim+1))
  {
      for(i=0; i<kdim; i++)
      {
          hvcoef[kk+i] = epcoef[kj+i] * evweight[ki]
			 + evcoef[kj+i] * epweight[ki];
      }

      hvcoef[kk+kdim] = evweight[ki];
  }


  /* Create C^1 cubic B-splines coefficients of the Hermite interpolant
     to the weights and control points. */

  s1379(scoef2,hvcoef,par_arr,inbcrv,kn1*(kdim+1),&qc,&kstat);
  if (kstat < 0) goto error;


  /* Create the surface */

  kind = 2;
  kcopy = 1;

  /* The surface is turned so that u is along the curves, v is across
     the curves in the lofting direction. */

  *rsurf = newSurf (kn1, qc->in, kord1, qc->ik, sknot1, qc->et, qc->ecoef,
		    kind, kdim, kcopy);
  if (*rsurf == SISL_NULL)
    goto err171;

  /* Release the curve object */

  if (qc != SISL_NULL)
    freeCurve (qc);


  /* Task done */

  goto out;


  /* Error in allocation. */

err101:
  *jstat = -101;
  s6err ("s1508", *jstat, kpos);
  goto out;

  /* Could not create surface. */

err171:
  *jstat = -171;
  s6err ("s1508", *jstat, kpos);
  goto out;

  /* Error in interpolation conditions. No. of curves < 2. */

err179:
  *jstat = -179;
  s6err ("s1508", *jstat, kpos);
  goto out;


  /* Error in lower level routine.  */

error:
  *jstat = kstat;
  s6err ("s1508", *jstat, kpos);
  goto out;

  /* Free allocated scratch  */

out:
  if (sknot1 != SISL_NULL)
    freearray (sknot1);
  if (scoef2 != SISL_NULL)
    freearray (scoef2);

  if (ratcurves != SISL_NULL)
  {
    for(i=0; i<inbcrv; i++)
    {
      if(ratcurves[i] != SISL_NULL) freeCurve(ratcurves[i]);
    }
    freearray(ratcurves);
  }
  if (epcoef != SISL_NULL) freearray(epcoef);
  if (evcoef != SISL_NULL) freearray(evcoef);
  if (epweight != SISL_NULL) freearray(epweight);
  if (evweight != SISL_NULL) freearray(evweight);
  if (epweight2 != SISL_NULL) freearray(epweight2);
  if (evweight2 != SISL_NULL) freearray(evweight2);
  if (hvcoef != SISL_NULL) freearray(hvcoef);

  return;
}
