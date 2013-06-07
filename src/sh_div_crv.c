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
 * $Id: sh_div_crv.c,v 1.2 2001-03-19 16:06:04 afr Exp $
 *
 */


#define SH_DIV_CRV

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
    sh_div_crv (SISLCurve * pc, int which_end, double aepsge, SISLCurve ** rcnew, int *jstat)
#else
void 
   sh_div_crv (pc, which_end, aepsge, rcnew, jstat)
     SISLCurve *pc;
     int which_end;
     double aepsge;
     SISLCurve **rcnew;
     int *jstat;
#endif
/*
********************************************************************
*
*********************************************************************
*
* PURPOSE     :To factorize a bezier curve J(x) over the interval [a,b]
*              into 
*              oldcurve = (x-a)/(b-a)*newcurve 
*              when which_end eq 0   (requires C0 = 0)
*                 and
*              oldcurve = (b-x)/(b-a)*newcurve 
*              when which_end eq 1   (requires Cn = 0).
*
*
*
*
* INPUT      : pc           - Oldcurve to factorize.
*              which_end     - Branch parameter for zero point.
*              aepsge       - Geometry tolerance.
*
*
*
* OUTPUT     : rcnew      -The new curve.
*              jstat     - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*              
*
* WRITTEN BY : Ulf J. Krystad, SI, 92-12.
* MODIFIED BY :
*
**********************************************************************/
{
  int kpos = 0;			/* Position of error.               */
  int ki,kj;                    /* Loop control                     */
  int kn,kk,kdim;               /* Attributes of inut curve         */			/* Position of error.               */
  double a,b;                   /* Bezier interval                  */
  double *et_new = SISL_NULL;        /* New knot array                   */
  double *ecoef_new = SISL_NULL;     /* New coefficient array            */
  SISLCurve *qc = SISL_NULL;		/* Pointer to new curve-object.     */


  /* Check that we have a curve. */
  if (!pc)
    goto err150;

  /* Minimum order allowed is 3. */
  if (pc->ik < 3)
     goto err151;
  
  /* The curve has to be of bezier type. */
  if (pc->in != pc->ik)
     goto err152;

  kn = pc->in;
  kk = pc->ik;
  a  = pc->et[kk-1];
  b  = pc->et[kn];
  kdim = pc->idim;

    /* Test if the corresponding coeficient is zero. */
/*  if (which_end == 0)
  {
     for (ki=0; ki < kdim; ki++)
	if (fabs(pc->ecoef[ki]) > aepsge)
	   goto err153;
  }
  else
  {
     for (ki=(kn-1)*kdim; ki < kn*kdim; ki++)
	if (fabs(pc->ecoef[ki]) > aepsge)
	   goto err153;
  }
  */
  
  /* create knot array. __________________________________________*/
  if ((et_new= newarray(kn+kk-2,DOUBLE)) == SISL_NULL) goto err101;

  for (ki=0; ki < kk-1; ki++)
  et_new[ki] = a;

    for (; ki < kn+kk-2; ki++)
  et_new[ki] = b;

  /* create coeficient array. _________________________________ */
  if ((ecoef_new= newarray(kdim*(kn-1),DOUBLE)) == SISL_NULL) goto err101;

  if (which_end)
     for (ki=0; ki < kn-1; ki++)
	for (kj=0; kj < kdim; kj++)
	   ecoef_new[ki*kdim +kj] = pc->ecoef[ki*kdim +kj]*(kn-1)/(kn-1-ki);
  else
     for (ki=0; ki < kn-1; ki++)
	for (kj=0; kj < kdim; kj++)
	   ecoef_new[ki*kdim +kj] = pc->ecoef[(ki+1)*kdim + kj]*(kn-1)/(ki+1);
  
  
  /* Create factor curve */
  if ((qc = newCurve (kn-1, kk-1, et_new, ecoef_new, pc->ikind, kdim, 2))
      == SISL_NULL) goto err101;

  *rcnew = qc;
  *jstat = 0;
  goto out;

/* ERROR EXITS ___________________________________________ */

/* Error. No input curve.  */
err150:
  *jstat = -150;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;


/* Error. order less than 3.  */
err151:
  *jstat = -151;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;

/* Error. Not a bezier curve.  */
err152:
  *jstat = -152;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;


/* Error in allocation.*/

err101:
  if (et_new) freearray(et_new);
  if (ecoef_new) freearray(ecoef_new);
  *jstat = -101;
  s6err ("sh_div_crv", *jstat, kpos);
  goto out;

out:
;
}
