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
 *
 *
 */


#define S1516

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1516(double ep[],double epar[],int im,int idim,
           double **ev,int *jstat)
#else
void s1516(ep,epar,im,idim,ev,jstat)
     double ep[];
     double epar[];
     int    im;
     int    idim;
     double **ev;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose:   To estimate the first derivative at each point in a sequence.
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          epar   - Parametrization array. The array should be increasing
*                   in value.
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
* Output:
*          ev     - Pointer to array containing the derivatives in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI 1993-10
*
*********************************************************************
*/
{
  int ki,kj;          /* Loop variables                              */
  int kk;             /* Polynomial order                            */
  int kpos=0;         /* Position of error                           */
  int kstat=0;        /* Status variable                             */
  double *gpar;
  int kcnsta;
  int kcnend;
  int iopen;
  int iorder;
  int ileft;
  double *ntype;
  SISLCurve *qc;
  int knbpar;
  double *evtemp;
  double cendpar;
  double *eder;




  /* Check input */

  if (idim < 1 || im < 2) goto err102;


  /* Allocate array for derivatives */

  evtemp    = newarray(idim*im,DOUBLE);
  if (evtemp == SISL_NULL) goto err101;

  ntype    = newarray(im,DOUBLE);
  if (ntype == SISL_NULL) goto err101;

  for(ki=0; ki<im; ki++)
  {
      ntype[ki] = 1.0;
  }

  eder    = newarray(2 * idim,DOUBLE);
  if (eder == SISL_NULL) goto err101;



  kcnsta = 1;
  kcnend = 1;
  iopen = 1;
  iorder = 4;

  s1358(ep, im, idim, ntype, epar, kcnsta, kcnend, iopen, iorder,
        epar[0],&cendpar, &qc, &gpar, &knbpar, &kstat);
     if(kstat < 0) goto error;

  for(ki=0,kk=0; ki<im; ki++,kk+=idim)
  {
      s1221(qc,1,epar[ki],&ileft,eder,&kstat);
      if(kstat < 0) goto error;

      for(kj=0; kj<idim; kj++)
      {
          evtemp[kk+kj] = eder[idim+kj];
      }
  }


  /* Calculation completed */

  /* Set result. */

  (*ev) = evtemp;

  *jstat = 0;
  goto out;



  /* Error in space allocation */

 err101: *jstat = -101;
  s6err("s1516",*jstat,kpos);
  goto out;


  /* Error in input. */

 err102: *jstat = -102;
  s6err("s1516",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */

 error:  *jstat =kstat;
  s6err("s1516",*jstat,kpos);
  goto out;

 out:
  if (ntype != SISL_NULL) freearray(ntype);
  if (gpar != SISL_NULL) freearray(gpar);
  if (eder != SISL_NULL) freearray(eder);

  return;
}
