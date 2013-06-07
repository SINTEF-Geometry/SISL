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
 * $Id: s1355.c,v 1.3 2001-03-19 15:58:47 afr Exp $
 *
 */

#define S1355

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1355(SISLCurve *pc,double eeps[],double **epar,
	   int *im, int *jstat)
#else
void s1355(pc,eeps,epar,im,jstat)
     SISLCurve  *pc;
     double eeps[];
     double **epar;
     int *im;
     int    *jstat;
#endif
/*
***********************************************************
*
* Purpose: To determine a set of sampling points on a given spline
*          curve such that the error in the cubic hermite spline
*          interpolant will be smaller than a given tolerance.
*
* Input :
*         Pc     - Pointer to spline curve.
*         Eeps   - Array (length kdim) giving the desired accuracy of
*                  the hermite interpolant in each component.
*
* Output:
*         Epar   - Array of dimension im containing the parametervalues
*                  of the desired sampling points.
*         Im     - Number of sampling points.
*         Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning:
*
* Method:
*
* Calls: s1707, s1720, s6err
*
* Written by: C.R. Birkeland, Si, April 1993.
* Changed by: Per OEyvind, SINTEF, 1994-11.
*             Added input check for invalid Eeps values.
**********************************************************
*/
{
  int i,j,k;               /* Loop control variables                       */
  int number = 0;
  int idim = pc->idim;     /* Space-dimension                              */
  int ik, in;              /* Order and number of vertices                 */
  int stat = 0;            /* Error control variables                      */
  int kpos = 0;
  int multi = 0;           /* Used to study multiplicity of original curve */
  int parlen;
  int help;
  double end, start;       /* Used to determine sampling point parameters  */
  double maxstep;
  double *coeff;           /* Only used as pointer to curve coefficients   */
  double *et = SISL_NULL;       /* Pointer to knotvector                        */
  double *par = SISL_NULL;      /* Storage of SISL-curve parameters             */
  double *maxtab = SISL_NULL;
  double *neweps = SISL_NULL;
  SISLCurve *oc=SISL_NULL;      /* Local SISL-curve                             */

  /* Check input-curve. */

  if (!pc) goto err150;
  s1707(pc, &stat);
  if (stat<0) goto error;

  /* Check other input */

  if (pc->in < pc->ik || pc->ik < 1 || idim < 1) goto err103;

  /* Find 4'th derivative of spline pc used in error estimate */

  s1720( pc, 4, &oc, &stat);
  if (stat<0) goto error;

  /* Initialize */

  ik = oc->ik;
  in = oc->in;
  et = oc->et;
  coeff = oc->ecoef;

  /* Allocate matrices */

  parlen = ik*in;
  par    = newarray(parlen, DOUBLE);
  maxtab = newarray(idim, DOUBLE);
  neweps = newarray(idim, DOUBLE);
  if (maxtab == SISL_NULL || par == SISL_NULL || neweps == SISL_NULL) goto err101;

  /* Do calculations here which may be done outside main loop */

  for(i=0; i<idim; i++)
    neweps[i] = eeps[i] * 384 / 5. ;

  /* Main loop */

  number = 0;
  multi = 0;
  for(help=0; help<in-ik+1; help++)
    {
      start = et[ik+help-1];
      end = et[ik+help];

      /* Find largest coefficient in each direction
       * to determine largest possible 4'th derivative
       * in interval et[start], et[end]                 */

      /* Set maxtab to zero */

      for(i=0; i<idim; i++)
	maxtab[i] = 0.;
      for(i=help*idim, k=0; k<ik-1; i+=idim, k++)
	for(j=0; j<idim; j++)
	  maxtab[j] = MAX( maxtab[j], fabs(coeff[i+j]) );

      /* Find maximum interval between two parametervalues,
       * in interval [ et(start), et(end) ].
       * Minimum step (Set a little bit larger than 'end - start' */

      maxstep = 1.1*end-start;
      for(i=0; i<idim; i++)
	{
	  if( maxtab[i] > 0. )
	    maxstep = MIN( neweps[i]/maxtab[i], maxstep );
	}
      maxstep = pow(maxstep, 0.25);

      /* Validate maxstep (and thus eeps) */

      if (DEQUAL(maxstep, 0.0))
        goto err103;

      /* Max. distance between parametervalues in interval
       * [start, end] must not exceed maxstep.
       * Generate the necessary parametervalues            */

      while(start<end)
	{
	  par[number] = start;
	  start += maxstep;
	  number++;

	  /* Make sure array for parametervalues is large enough */

	  if (number+1>=parlen)
	    {
	      parlen = 2*parlen;
	      if ((par = increasearray(par, parlen, double)) == SISL_NULL)
		goto err101;
            }
	}

      /* Test for multiple knots in original knotvector */

      if(end == start)
	{
	  /* Check if knotmultiplicity is  'pc->ik - 1' */

	  if(multi == pc->ik - 3)
	    {
	      par[number] = start;
	      number++;
	      multi=0;
	    }
	  else
	    multi++;
	}
      else
	multi = 0;
    }

  par[number] = end;
  if ((par = increasearray(par,number+1,DOUBLE)) == SISL_NULL) goto err101;
  *epar = par;
  *im = number+1;

  /* Success ! */

  *jstat = 0;
  goto out;

  /* Allocation error. */

 err101:
   *jstat = -101;
   s6err("s1355",*jstat,kpos);
   goto out;

  /* Error in input */

  err103:
    *jstat = -103;
    s6err("s1355",*jstat,kpos);
    goto out;

  /* Empty curve. */

  err150:
    *jstat = -150;
    s6err("s1355",*jstat,kpos);
    goto out;

  /* Error in lower level routine. */

  error:
    *jstat = stat;
    s6err("s1355",*jstat,kpos);
    goto out;

  /* Exit */

  out:
    if( maxtab != SISL_NULL) freearray(maxtab);
    if( neweps != SISL_NULL) freearray(neweps);
    if( oc != SISL_NULL) freeCurve(oc);
    return;
}
