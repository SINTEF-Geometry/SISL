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
 * $Id: s1750.c,v 1.2 2001-03-19 15:58:53 afr Exp $
 *
 */


#define S1750

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1750(SISLCurve *pc,int ikh,SISLCurve **rc,int *jstat)
#else
void s1750(pc,ikh,rc,jstat)
     SISLCurve *pc;
     int   ikh;
     SISLCurve **rc;
     int   *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a B-spline curve using a higher order
*              B-spline basis.
*             
* INPUT      : pc     - The input B-spline curve.
*              ikh    - Order of the new urve.
*
* OUTPUT     : rc     - Pointer to the higher order curve
*              jstat  - status messages :
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : The order of the curve is elevated one order at the time
*              using the algorithm of COHEN, LYCHE and SCHUMAKER.
*
* REFERENCES :
*
*
* CALLS      : s1754,s1753,s1755,s6err.
*
*
* WRITTEN BY : 	Qyvind Hjelle, SI, Oslo, Norway. 10. Nov 1988
* REWRITTEN BY:	Christophe R. Birkeland, SI, 1991-07
* REVISED BY :  Johannes Kaasa, SI, May 1992 (Introduced NURBS)
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/
{
  int ki, kn, kk;		/* Loop control parameters			*/
  int kordr;
  int inhrem;			/* Used to store inh, for later use in last
				 * call to s1753				*/
  int kpos = 0;			/* Error position indicator			*/
  int kstat = 0;		/* Status variable */
  double *kcc = SISL_NULL;
  double *kcw = SISL_NULL;		/* Arrays for internal use only			*/
  double *orknot = SISL_NULL;	/* Used to store 'original' knot vector		*/
  double *xtknot = SISL_NULL;	/* Used to store extended knot vector		*/
  double *pointer = SISL_NULL;
  double *orcoef = SISL_NULL;	/* Used to store 'original' coefficient matrix	*/
  double *et = SISL_NULL;		/* Original knot vector				*/
  double *ebcoef = SISL_NULL;	/* Vertices of original curve			*/
  int in;			/* Number of vertices of original curve		*/
  int ik;			/* Order of original curve			*/
  int idim;			/* Dimension of the space where the curve lie	*/
  int kdim;                     /* Potential rational dimension.                */
  int kind;                     /* Kind of curve, 2 and 4 are rationals.        */
  double *iknt = SISL_NULL;		/* New knot vector				*/
  double *icoef = SISL_NULL;		/* Coefficients of new curve			*/
  int inh;			/* Number of vertices produced			*/

  *jstat = 0;

  /* Initialization of variables. */

  kind = pc->ikind;
  idim = pc->idim;
  et = pc->et;
  if (kind == 2 || kind == 4)
    {
      ebcoef = pc->rcoef;
      kdim = idim + 1;
    }
  else
    {
      ebcoef = pc->ecoef;
      kdim = idim;
    }
  in = pc->in;
  ik = pc->ik;

  /* Test if legal input. */

  if ((ik < 1) || (ikh < ik) || (in <ik)) goto err112;

  /* If ikh=ik, copy input curve to output variables. */

  if (ikh == ik)
    {
      *rc = newCurve (in, ik, et, ebcoef, pc->ikind, idim, 1);
      if (*rc == SISL_NULL) goto err171;

      /* If the input curve is periodic, the output curve is periodic. */
      (*rc)->cuopen = pc->cuopen;
      goto out;
    }

  /* Find size of knot vector and vertex vector,
     and find knot vector expressed in order ikh. */

  s1754 (et, in, ik, ikh, &iknt, &inh, &kstat);
  if (kstat < 0) goto error;

  /* Allocate coefficients array for raised curve. */

  if((icoef = newarray (inh * kdim, DOUBLE)) == SISL_NULL) goto err101;

  /* Allocate arrays for internal use. */

  if((kcc = newarray (kdim * ikh, DOUBLE)) == SISL_NULL) goto err101;
  if((kcw = newarray (kdim * ikh, DOUBLE)) == SISL_NULL) goto err101;

  /* Find vertices if  ikh = ik+1 */

  if (ikh == ik + 1)
    {
      s1753 (et, ebcoef, in, ik, kdim, iknt, icoef, inh, kcc, kcw, &kstat);
      if (kstat < 0) goto error;

      *rc = newCurve (inh, ikh, iknt, icoef, pc->ikind, idim, 2);
      if (*rc == SISL_NULL) goto err171;

      /* If the input curve is periodic, the output curve is periodic. */
      (*rc)->cuopen = pc->cuopen;
 
      goto out;
    }

  /* Allocate arrays to store knot vector for use in s1755. */

  orknot = newarray ((in +ik) *(ikh - ik + 1), DOUBLE);
  if (orknot == SISL_NULL) goto err101;
  xtknot = newarray ((in +ik) *(ikh - ik + 1), DOUBLE);
  if (xtknot == SISL_NULL) goto err101;

  /* Allocate array to store vertices. */

  orcoef = newarray (inh * kdim, DOUBLE);
  if (orcoef == SISL_NULL) goto err101;

  /* Initialize orknot and orcoef. */

  for (ki = 0; ki < (in +ik); ki++)
    orknot[ki] = et[ki];

  for (ki = 0; ki < (kdim * in); ki++)
    orcoef[ki] = ebcoef[ki];


  /* MAIN LOOP. Do the order raisings. */

  inhrem = inh;
  kn = in;
  kk = ik;
  for (kordr = ik + 1; kordr < ikh; kordr++)
    {
      /* Produce raised knots. */

      s1755 (orknot, kn, kk, xtknot, &inh, &kstat);
      if (kstat < 0) goto error;

      /* Produce raised vertices. */

      s1753 (orknot, orcoef, kn, kk, kdim, xtknot, icoef,
	     inh, kcc, kcw, &kstat);
      if (kstat < 0) goto error;


      if ((kordr + 1) < ikh)
	{
	  pointer = orknot;
	  orknot = xtknot;
	  xtknot = pointer;
	}
      kk = kordr;
      kn = inh;
      pointer = orcoef;
      orcoef = icoef;
      icoef = pointer;
    }

  inh = inhrem;
  s1753 (xtknot, orcoef, kn, kk, kdim, iknt, icoef, inh, kcc, kcw, &kstat);
  if (kstat < 0) goto error;

  /* OK.
   * Create new curve */

  *rc = newCurve (inh, ikh, iknt, icoef, pc->ikind, idim, 2);
  if (*rc == SISL_NULL) goto err171;

  /* If the input curve is periodic, the output curve is periodic. */
  (*rc)->cuopen = pc->cuopen;

  goto out;


  /* Error in array allocation */

  err101:
    *jstat = -101;
    s6err ("s1750", *jstat, kpos);
    goto out;

  /* Could not create curve. */

  err171:
    *jstat = -171;
    s6err ("s1750", *jstat, kpos);
    goto out;

  /* Error in description of B-spline */

  err112:
    *jstat = -112;
    s6err ("s1750", *jstat, kpos);
    goto out;

  /* Error in lower level routine */

  error:
    *jstat = kstat;
    s6err ("s1750", *jstat, kpos);
    goto out;

  out:
    if (kcc != SISL_NULL)    freearray (kcc);
    if (kcw != SISL_NULL)    freearray (kcw);
    if (orknot != SISL_NULL) freearray (orknot);
    if (xtknot != SISL_NULL) freearray (xtknot);
    if (orcoef != SISL_NULL) freearray (orcoef);
    return;
}
