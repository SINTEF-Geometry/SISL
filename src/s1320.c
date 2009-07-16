/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s1320.c,v 1.2 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1320

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1320 (SISLSurf * psurf, double earray[], int inarr,
       int ratflag, SISLSurf ** rsurf, int *jstat)
#else
void
s1320 (psurf, earray, inarr, ratflag, rsurf, jstat)
     SISLSurf *psurf;
     double earray[];
     int inarr;
     int ratflag;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To put a surface description into the implicit
*              second order surface described by the input array.
*
* INPUT      : psurf  - Pointer to input surface.
*              earray - The description of the input array
*                       dimension (psurf->idim+1)^2*inarr
*              inarr  - Number of parallel matrices in earray.
*                       inarr should be less or equal to 3.
*              ratflag - If psurf is nonrational it is ignored.
*                        Otherwise:
*                        If ratflag = 0 rsurf is the nonrational numerator
*                        If ratflag = 1 rsurf is a full rational surface
*
* OUTPUT     : rsurf  - The resulting surface
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
* METHOD     : Dependent on the type of object we make:
*
*        F(S,T) = (P(s,t),1)  EARRAY (P(s,t),1)
*
*     by sampling enough points to use interpolation for reproduction.
*
* REFERENCES :
*
* CALLS      : s1896,s6err.
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway.
* REVISED BY : Mike Floater, 91-01, SI, Oslo, Norway for rational surfaces.
* REVISED BY : Mike Floater, SI, Oslo 11/9/91 -- ratflag.
* REVISED BY : Michael Floater, SI, June 92. The rational stuff
*              was messed up after translation of
*              the fortran part to c. But it works now.
* REVISED BY : Atgeirr F Rasmussen, Sintef, October 2000. Fixed the
*              rational stuff, which assumed dimension was 3.
*
*********************************************************************
*/
{
  int kpos = 0;
  int kstat = 0;
  SISLSurf *ssurf = SISL_NULL;	/* Temperary SISL-surface. */
  int kdim;			/* Number of dimesions in psurf                     */
  int kdimp1;			/* Dimension of  earray should be kdim+1            */
  int lder[3];			/* Derivative indicator array                       */
  double *scoef = SISL_NULL;		/* Vertices of psurf (scaled in the rational case)  */
  double *rscoef = SISL_NULL;	/* pointer to vertices in the rational case         */
  int ikind;			/* kind of surface                                  */
  double wmin, wmax;		/* min. and max. weight values for rational surface */
  double scale;			/* factor used for scaling rational weights         */
  int i;			/* loop variable                                    */
  double *sarray = SISL_NULL;	/* Array for calculating denominator if used      */
  int knarr;			/* Number of parallel arrays to use.   */
  int nkind;			/* Kind of output surface (rsurf).    */
  SISLSurf *jsurf = SISL_NULL;       /* Temporary SISLSurf. */

  *jstat = 0;


  /* Make local pointers. */

  kdim = psurf->idim;
  ikind = psurf->ikind;

  /* Set dimension of kdimp1.  */

  kdimp1 = kdim + 1;


  /* Test input. */

  if (kdim < 1)
    goto err102;
  if (inarr < 1 || 3 < inarr)
    goto err172;


  /* rational surfaces is a special case. */

  if (ikind == 2 || ikind == 4)
    {
      kdim++;
      /* scale the coeffs so that min. weight * max. weight = 1. */

      rscoef = psurf->rcoef;
      wmin = rscoef[kdim-1];
      wmax = rscoef[kdim-1];

      for (i = kdim-1; i < psurf->in1 * psurf->in2 * kdim; i += kdim)
	{
	  if (rscoef[i] < wmin)
	    wmin = rscoef[i];
	  if (rscoef[i] > wmax)
	    wmax = rscoef[i];
	}

      scale = (double) 1.0 / sqrt (wmin * wmax);
      scoef = newarray (psurf->in1 * psurf->in2 * kdim, DOUBLE);
      if (scoef == SISL_NULL)
	goto err101;

      for (i = 0; i < psurf->in1 * psurf->in2 * kdim; i++)
	{
	  scoef[i] = rscoef[i] * scale;
	}
    }
  else
    {
      scoef = psurf->ecoef;
    }

  ssurf = newSurf (psurf->in1, psurf->in2, psurf->ik1, psurf->ik2,
		   psurf->et1, psurf->et2, scoef, 1, kdim, 1);
  if (ssurf == SISL_NULL)
    goto err171;

  if ((ikind == 2 || ikind == 4) && ratflag == 1)
    {
      /* Output surface will also be rational. */

      nkind = 2;

      /* Add an extra parallel array to pick up the weights
	 of the subsequent homogeneous vertices of rsurf. */

      knarr = inarr + 1;

      sarray = new0array (kdimp1 * kdimp1 * knarr, DOUBLE);
      if (sarray == SISL_NULL)
	goto err101;

      memcopy (sarray, earray, kdimp1 * kdimp1 * inarr, double);

      sarray[kdimp1 * kdimp1 * knarr - 1] = (double) 1.0;
    }
  else
    {
      nkind = 1;
      knarr = inarr;
      sarray = earray;
    }

  lder[0] = 0;
  lder[1] = 0;
  lder[2] = 0;

  /* Put surface into implicit surface */

  s1896 (ssurf, sarray, kdimp1, knarr, lder, lder, lder, lder, &jsurf, &kstat);
  if (kstat < 0)
    goto error;

  if ((ikind == 2 || ikind == 4) && ratflag == 1)
    {
      /* Output from s1896 is a dim+1 non-rational surface jsurf. */
      /* Convert homogeneous jsurf to rational rsurf. */

      *rsurf = newSurf(jsurf->in1,jsurf->in2,
                        jsurf->ik1,jsurf->ik2,
                        jsurf->et1,jsurf->et2,
                        jsurf->ecoef,
                        2,jsurf->idim-1,1);
      freeSurf(jsurf);
    }
  else
    {
      *rsurf = jsurf;
    }

  if (ikind == 2 || ikind == 4)
    {
      if (scoef)
	freearray (scoef);
      if (ratflag)
	freearray (sarray);
    }

  /* Ok. */

  goto out;


  /* Error in lower level function */

error:
  *jstat = kstat;
  s6err ("s1320", *jstat, kpos);
  goto out;

  /* allocation problems. */
err101:
  *jstat = -101;
  s6err ("s1320", *jstat, kpos);
  goto out;

  /* Dimension less than 1    */

err102:
  *jstat = -102;
  s6err ("s1320", *jstat, kpos);
  goto out;

  /* Could not create surface. */

err171:
  *jstat = -171;
  s6err ("s1320", *jstat, kpos);
  goto out;

  /* Dimension inarr not equal to 1,2 or 3 */

err172:
  *jstat = -172;
  s6err ("s1320", *jstat, kpos);
  goto out;

out:
  if (ssurf)
    freeSurf (ssurf);
  return;
}
