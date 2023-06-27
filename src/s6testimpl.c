/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF Digital,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF Digital, Department of Mathematics and Cybernetics,                         
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
 * written agreement between you and SINTEF Digital. 
 */

#include "sisl-copyright.h"


#define S6TESTIMPL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
static int s6testimpl_classify(SISLSurf *ps1, SISLSurf *ps2, SISLIntpt *vintpt[],
			       int inmbpt, int ndir[], int *jstat);
#else
static int s6testimpl_classify();
#endif


#if defined(SISLNEEDPROTOTYPES)
void
s6testimpl(SISLSurf *ps1, SISLSurf *ps2, int first, SISLIntpt *vintpt[],
	   int inmbpt, double aepsge, int *jstat)
#else
void
  s6testimpl(ps1, ps2, first, vintpt, inmbpt, aepsge, jstat)
   SISLSurf *ps1;
   SISLSurf *ps2;
   int first;
   SISLIntpt *vintpt[];
   int inmbpt;
   double aepsge;
   int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Try to intercept intersection between two B-spline
*              surfaces by finding a computing an implicit surface
*              approximating the first surface and inserting the
*              second surface into the implicit expression. Box test
*              and simple case test is performed.
*
*
* INPUT      : ps1      - First surface.
*              ps2      - Second surface.
*              first    - Indicates if first surface is to be implicitized
*              vintpt   - Intersection points on edges
*              inmbpt   - Number of intersection points on edges
*              aepsge   - Geometry tolerance. 
*
*
* OUTPUT     : jstat    - status messages
*                                = 3   : Simple case encountered
*                                = 1   : Intersection still possible.
*                                = 0   : Ok. No intersection is possible.
*                                < 0   : Error.
*
*
* METHOD     :
*
* CALLS      : 
*
* REFERENCES :
*
* WRITTEN BY : Vibeke Skytt, SINTEF, 2023-03.
*
*********************************************************************
*/
{
   int kstat = 0;   /* Local status variable.  */
   int kdim = ps1->idim;  /* Dimension of space. */
   double tang;     /* Angle between normal cone axes */
   int deg = 3;     /* Implicit degree                */
   double tepsge;   /* Local tolerance.        */
   SISLSurf *qs1 = SISL_NULL; /* 1D surface.                     */
   SISLSurf *qs2 = SISL_NULL; /* 1D surface.                     */
   int kk1, kk2;    /* Orders of first 1D surface.       */
   int kk3, kk4;    /* Orders of second 1D surface.      */
   int kn1, kn2;    /* Number of vertices of 1. 1D surface. */
   int kn3, kn4;    /* Number of vertices of 2. 1D surface. */
   double *st1 = NULL;   /* Knot vector of 1. 1D surface in 1. par. dir. */
   double *st2 = NULL;   /* Knot vector of 1. 1D surface in 2. par. dir. */
   double *st3 = NULL;   /* Knot vector of 2. 1D surface in 1. par. dir. */
   double *st4 = NULL;   /* Knot vector of 2. 1D surface in 2. par. dir. */
   double *sc = NULL;    /* Matrix of equation system for unknown implicit
			    coefficients */
   double *sc1 = NULL;   /* Coefficients of 1D surface.     */
   double *sc2 = NULL;   /* Coefficients of 1D surface.     */
   double *scimp = NULL; /* Implicit coefficents.           */
   int numscimp = 0;
   double grad;          /* Estimate of medium gradient.               */
   SISLSurf *ps3 = (first) ? ps1 : ps2;
   SISLSurf *ps4 = (first) ? ps2 : ps1;
   int ldir[10];
   int conn = 0;
   int ki;

   /* Check if the surfaces are of Bezier type */
   if (ps1->in1 > ps1->ik1 || ps1->in2 > ps1->ik2 ||
       ps2->in1 > ps2->ik1 || ps2->in2 > ps2->ik2)
     {
       *jstat = 1;
       goto out;
     }

   /* Check for normal cones in approximately the same area */
   tang = s6ang(ps1->pdir->ecoef,ps2->pdir->ecoef,kdim);
   if (tang > 0.25*PI)
   {
      *jstat = 1;
      goto out; 
   } 
   
   *jstat = 1;  /* Initiate to possibility of intersection.  */

   /* Make matrix if implicitization equation */
   s1326(ps3, deg, NULL, 1, &st1, &st2, &sc, &kk1, &kk2, &kn1, &kn2,
	 &numscimp, &kstat);
   if (kstat < 0)
     goto error;
   
   if (st1) freearray(st1);
   st1 = NULL;
   if (st2) freearray(st2);
   st2 = NULL;

   /* Solve underdetermined equation system to find implicit coefficents. */
   if ((scimp = newarray(numscimp, DOUBLE)) == NULL)
     goto err101;
   s1339(ps3, deg, sc, numscimp, kn1*kn2, aepsge, scimp, &grad, &kstat);
   if (kstat < 0)
     goto error;
   if (kstat == 1)
     goto warn1;  /* No implicit approximation computed. */
   if (sc != NULL) freearray(sc);
   sc = NULL;

   /* Local tolerance */
   tepsge = fabs(grad)*aepsge;

   /* Put the description of the surfaces into the implicit equation */
    s1326(ps3, deg, scimp, 1, &st1, &st2, &sc1, &kk1, &kk2, &kn1, &kn2,
	 &numscimp, &kstat);
   if (kstat < 0)
     goto error;
   
   if ((qs1 = newSurf(kn1, kn2, kk1, kk2, st1, st2, sc1, 1, 1, 2)) == NULL)
      goto err101;
   
   s1326(ps4, deg, scimp, 1, &st3, &st4, &sc2, &kk3, &kk4, &kn3, &kn4,
	 &numscimp, &kstat);
   if (kstat < 0)
     goto error;
   
   if ((qs2 = newSurf(kn3, kn4, kk3, kk4, st3, st4, sc2, 1, 1, 2)) == NULL)
      goto err101;
	 
   /* Make box of first 1D surface. */
   sh1992su(qs1,0,tepsge,&kstat);
   if (kstat < 0)
     goto error;
   sh1992su(qs1,2,tepsge,&kstat);
   if (kstat < 0)
     goto error;
   
   /* Make box of second 1D surface. */
   sh1992su(qs2,0,tepsge,&kstat);
   if (kstat < 0)
     goto error;
   sh1992su(qs2,2,tepsge,&kstat);
   if (kstat < 0)
     goto error;
   
   /* Check if the boxes overlap.  */
   if (qs1->pbox->e2min[2][0] > qs2->pbox->e2max[2][0] ||
       qs1->pbox->e2max[2][0] < qs2->pbox->e2min[2][0])
     {
       if (inmbpt > 0 &&
	   (qs1->pbox->e2min[0][0] <= qs2->pbox->e2max[0][0] &&
	    qs1->pbox->e2max[0][0] >= qs2->pbox->e2min[0][0]))
	 {
	   /* Extra test on point configuration */
	   if (inmbpt <= 10)
	     {
	       /* Classify points */
	       conn = s6testimpl_classify(ps1, ps2, vintpt, inmbpt, ldir, &kstat);
	       if (kstat < 0)
		 goto error;
	       for (ki=0; ki<inmbpt; ++ki)
		 if (ldir[ki] == -1 || ldir[ki] == 1)
		   break;   /* Intersection curve tangens points in or out */
	       if (conn || ki == inmbpt)
		 {
		   *jstat = 0;  /* No more intersections are possible */
		   goto out;
		 }
	       
	     }
	 }
       else
	 {
	   *jstat = 0;  /* No intersection possible */
	   goto out;
	 }
     }

   /* Perform simple case test. */
   if (ps1->pdir->igtpi != 0 || ps2->pdir->igtpi != 0)
     goto warn1;
   if (max(fabs(qs1->pbox->e2min[0][0]), fabs(qs1->pbox->e2max[0][0])) > tepsge)
     goto warn1;   /* Check with non-expanded box. 
		      The implicit approximation is not accurate enough
		      for a reliable simple case test */
   sh1994(qs2, tepsge, &kstat);
   *jstat = (kstat == 1) ? 3 : 1;
   goto out;
   
   warn1 : *jstat = 1; /* No implicit approximation computed or 
			  complex surface. */
   goto out;           /* Possibility of intersection.        */
   
   error : *jstat = kstat;
   goto out;
   
   err101 : *jstat = -101;         /* Error in scratch allocation. */
   s6err("s6testimpl",*jstat,0);
   goto out;
   
   out:
      if (qs1) freeSurf(qs1);
      if (qs2) freeSurf(qs2);
      if (sc) freearray(sc);
      if (scimp) freearray(scimp);
      
      return;
}


#if defined(SISLNEEDPROTOTYPES)
static int
s6testimpl_classify(SISLSurf *ps1, SISLSurf *ps2, SISLIntpt *vintpt[],
		    int inmbpt, int ndir[], int *jstat)
#else
static int
  s6testimpl_classify(ps1, ps2, vintpt, inmbpt, ndir, jstat)
     SISLSurf *ps1;
     SISLSurf *ps2;
     SISLIntpt *vintpt[];
     int inmbpt;
     int ndir[];
     int *jstat;
#endif
 /*
 *********************************************************************
 *
 *********************************************************************
 *
 * PURPOSE    : To classify intersection points on edges on two surfaces
 *
 *
 *
 * INPUT      : ps1      - First surface in intersection.
 *              ps2      - Second surface in intersection.
 *              vintpt   - Intersection points on edges
 *              inmbpt   - Number of intersection points
 *
 * OUTPUT     : retur    - Whether (1) or not (0) the edge points are connected
 *              ndir     - Classification of intersection curve in each point
 *                               * 0 - The intersect.curve is parallel to one
				 *      parameter direction.
				 *  1 - The intersect.curve has direction into the
				 *      domain.
				 * -1 - The intersect.curve has direction out of the
				 *      domain.
				 *  2 - The point is singulear.
				 * 10 - The intersect.curve touch one corner of the
				 *      domain.
 *              jstat    - status messages
 *                           = 0     : OK.
 *                           < 0     : error
 *
 *
 * METHOD     :
 *
 *
 * REFERENCES :
 *
 *
 * WRITTEN BY : Vibeke Skytt, SINTEF, 23-03
 *
 *********************************************************************
 */
{
  int kstat = 0;
  int ki, kj, kn, kv;
  int klist1, klist2;
  int klfs1 = 0, klft1 = 0, klfs2 = 0, klft2 = 0;
  int conn = 0;
  int dim = ps1->idim;
  double sstart[4], send[4];
  double tang;
  int kdir;
  double *sval1 = SISL_NULL;
  double *sval2, *snorm1, *snorm2, *stang, *sdec1, *sdec2;
  double *stmp, *stmp2;
  int ix;
  double tolpar = (double)0.001;

  *jstat = 0;

  sstart[0] = ps1->et1[ps1->ik1-1];
  send[0] = ps1->et1[ps1->in1];
  sstart[1] = ps1->et2[ps1->ik2-1];
  send[1] = ps1->et2[ps1->in2];
  sstart[2] = ps2->et1[ps2->ik1-1];
  send[2] = ps2->et1[ps2->in1];
  sstart[3] = ps2->et2[ps2->ik2-1];
  send[3] = ps2->et2[ps2->in2];
  
  /* Check if the intersection points are connected */
  for (ki=0; ki<inmbpt; ++ki)
    {
      for (kj=ki+1; kj<inmbpt; ++kj)
	{
	  sh6getlist(vintpt[ki], vintpt[kj], &klist1, &klist2, &kstat);
	  if (kstat < 0)
	    goto error;
	  if (kstat == 0)
	    break;   /* Connection */
	}
      if (kj == inmbpt)
	break;  /* No connection is found */
    }
  if (ki == inmbpt)
    conn = 1;

  /* Classify points */
  if ((sval1 = newarray (11*dim, double)) == SISL_NULL)
	goto err101;
  sval2 = sval1 + 3*dim;
  snorm1 = sval2 + 3*dim;
  snorm2 = snorm1 + dim;
  stang = snorm2 + dim;
  sdec1 = stang + dim;
  sdec2 = sdec1 + dim;

  for (ki=0; ki<inmbpt; ++ki)
    {
      ndir[ki] = -20;  /* Direction not set */
      s1421(ps1, 1, vintpt[ki]->epar, &klfs1, &klft1, sval1, snorm1, &kstat);
      if (kstat < 0)
	goto error;
      
       s1421(ps2, 1, vintpt[ki]->epar+2, &klfs2, &klft2, sval2, snorm2, &kstat);
      if (kstat < 0)
	goto error;

      tang = s6ang (snorm1, snorm2, dim);
      if (tang < REL_PAR_RES)
	{
	  ndir[ki] = 2;  /* Singularity */
	  continue;
	}

      s6crss (snorm1, snorm2, stang);

      s6decomp (stang, sdec1, sval1+dim, sval1+2*dim, snorm1, &kstat);
      if (kstat < 0)
	goto error;
      else if (kstat > 0)
	{
	  ndir[ki] = 2;
	  continue;
	}

      s6decomp (stang, sdec2, sval2+dim, sval2+2*dim, snorm2, &kstat);
      if (kstat < 0)
	goto error;
      else if (kstat > 0)
	{
	  ndir[ki] = 2;
	  continue;
	}

      for (int kj=0; kj<4; ++kj)
	{
	  kdir = -20;
	  stmp = (kj < 2) ? sval1 : sval2;
	  stmp2 = (kj < 2) ? sdec1 : sdec2;
	  ix = (kj == 0 || kj == 2) ? 2 : 1;
	  if (DEQUAL(vintpt[ki]->epar[kj], sstart[kj])) 
	    {
	      tang = s6ang(stang, stmp+ix*dim, dim);
	      kdir = (stmp2[2-ix] < DZERO) ? 1 : -1;
	    }
	  else if (DEQUAL(vintpt[ki]->epar[kj], send[kj]))
	    {
	      tang = s6ang(stang, stmp+ix*dim, dim);
	      kdir = (stmp2[2-ix] < DZERO) ? -1 : 1;
	    }
	  if (kdir != -20 && tang < tolpar)
	    kdir = 0;
	  if (ndir[ki] == -20)
	    ndir[ki] = kdir;
	  else if (ndir[ki] != kdir)
	    {
	      if (ndir[ki] == 0)
		ndir[ki] = kdir;
	      else
		{
		  ndir[ki] = 10;
		  break;
		}
	    }
	}
    }

  goto out;
  
  /* Error in sub rutines.      */
 error:*jstat = kstat;
  goto out;

  /* Error in memory allocation.      */
 err101:*jstat = -101;
  goto out;

 out:
  if (sval1 != SISL_NULL)
    freearray(sval1);

  return conn;
}
