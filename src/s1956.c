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
 * $Id: s1956.c,v 1.3 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1956

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1956(SISLCurve *pc1,SISLCurve *pc2,SISLSurf **rsurf,int *jstat)
#else
void s1956(pc1,pc2,rsurf,jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     SISLSurf  **rsurf;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Find the surface describing the difference function
*              between two B-spline curves.
*
*
* INPUT      : pc1    - Pointer to first curve in difference function.
*              pc2    - Pointer to second curve in difference function.
*
*
*
* OUTPUT     : rsurf  - Pointer to difference surface.
*              jstat  - status messages
*                     = 2      : The curves have the same representation
*                                and opposite direction.
*                     = 1      : The curves have the same representation
*                                and the same direction.
*                     = 0      : ok
*                     < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s6dist - Distance between two points.
*              s6diff - Difference vector between two vectors.
*              newSurf - Create new surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
* REVISED BY : Michael Floater, SI, 91-09 for possible rational curves.
*          Also fixed a bug in the "check same representation" part.
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Sept. 1994. Fixed memory
*          leak from 'ratpc'.
*
*********************************************************************
*/
{
  int kstat = 0;        /* Local status variable.                        */
  int kpos = 0;         /* Position of error.                            */
  int i,ki,kj;          /* Counters.                                     */
  int kind1,kind2;      /* Kinds of the two curves.                       */
  int kdim;             /* Dimension of the space in which the curves lie.*/
  int kn1;              /* Number of vertices of first curve.            */
  int kk1;              /* Order of first curve.                         */
  int kn2;              /* Number of vertices of second curve.           */
  int kk2;              /* Order of second curve.                        */
  double tdist;         /* Distance between vertices.                    */
  double *scoef1;       /* Pointer to vertices of first curve.           */
  double *scoef2;       /* Pointer to vertices of second curve.          */
  double *ssurf = SISL_NULL; /* Pointer to the vertices of the difference
			   surface.                                      */
  double *s0;           /* Pointer used to traverse ssurf.               */
  double *s1;           /* Pointer used to traverse scoef1.              */
  double *s2;           /* Pointer used to traverse scoef2.              */
  SISLCurve *ratpc=SISL_NULL; /* Temporary SISLCurve for rational cases.      */
  int newdim;          /* = kdim+1 if pc1 or pc2 rational, else = kdim.  */
  double *s1weight;    /* Pointer to weights in pc1 if rational.         */
  double *s2weight;    /* Pointer to weights in pc2 if rational.         */
  double *s3;          /* Pointer to traverse scoef2 if rational.        */

  /* Test input. */

  if (pc1->idim != pc2->idim) goto err106;

  /* Describe curves in local parameters.  */

  kdim = pc1 -> idim;
  kn1 = pc1 -> in;
  kk1 = pc1 -> ik;
  kind1 = pc1 -> ikind;

  kn2 = pc2 -> in;
  kk2 = pc2 -> ik;
  kind2 = pc2 -> ikind;

  if((kind1 == 2 || kind1 == 4) || (kind2 == 2 || kind2 == 4))
  {
    /* At least one of the curves is rational. Convert the
       other curve to rational if necessary. */

    newdim=kdim+1;

    if(kind1 != 2 && kind1 != 4)
    {
	ratpc=s1521(pc1,&kstat);
	if(kstat < 0) goto err102;

	scoef1 = ratpc -> rcoef;
	scoef2 = pc2 -> rcoef;
    }
    else if(kind2 != 2 && kind2 != 4)
    {
	ratpc=s1521(pc2,&kstat);
	if(kstat < 0) goto err102;

	scoef1 = pc1 -> rcoef;
	scoef2 = ratpc -> rcoef;
    }
    else
    {
	scoef1 = pc1 -> rcoef;
	scoef2 = pc2 -> rcoef;
    }


    /* Allocate space for homogeneous coefficients of rational
       difference surface.  */

    if ((ssurf = newarray(kn1*kn2*newdim,double)) == SISL_NULL) goto err101;

    /* Make coefficients of difference function.  */

    for (s0=ssurf,s2=scoef2,kj=0; kj<kn2; s2+=(kdim+1),kj++)
    {
      for (s1=scoef1,ki=0; ki<kn1; s0++,s1++,ki++)
      {
	  s1weight=s1+kdim;
	  s2weight=s2+kdim;

	  for(i=0,s3=s2; i<kdim; i++,s0++,s1++,s3++)
	  {
	      (*s0) = (*s2weight) * (*s1) - (*s1weight) * (*s3);
	  }

	  (*s0) = (*s1weight) * (*s2weight);
      }
    }

    /* Create difference function.  */

    *rsurf = SISL_NULL;
    if ((*rsurf = newSurf(kn1,kn2,kk1,kk2,pc1->et,pc2->et,
			  ssurf,2,kdim,1)) == SISL_NULL) goto err101;
  }
  else
  {
      /* Both curves are non-rational. */

    newdim=kdim;

    scoef1 = pc1 -> ecoef;
    scoef2 = pc2 -> ecoef;

    /* Allocate space for coefficients of difference surface.  */

    if ((ssurf = newarray(kn1*kn2*kdim,double)) == SISL_NULL) goto err101;

    /* Make coefficients of difference function.  */

    for (s0=ssurf,s2=scoef2,kj=0; kj<kn2; s2+=kdim,kj++)
      for (s1=scoef1,ki=0; ki<kn1; s0+=kdim,s1+=kdim,ki++)
        s6diff(s1,s2,kdim,s0);

    /* Create difference function.  */

    *rsurf = SISL_NULL;
    if ((*rsurf = newSurf(kn1,kn2,kk1,kk2,pc1->et,pc2->et,
			  ssurf,1,kdim,1)) == SISL_NULL) goto err101;
  }



  /* Test if the curves have the same representation and same direction. */

  kstat = 1;
  if (kn1 != kn2 || kk1 != kk2) kstat = 0;
  tdist = s6dist(scoef1,scoef2,newdim);
  if(DNEQUAL(tdist,(double)0.0))
  {
      kstat = 0;
  }
  else
  {
      for (s1=scoef1+newdim,s2=scoef2+newdim,ki=1; ki<kn1 && kstat>0;
           s1+=newdim,s2+=newdim,ki++)
      {
        if (DNEQUAL(s6dist(s1,s2,newdim),tdist))
	{
	    kstat = 0;
	    break;
	}
      }
  }

  if (kstat == 0)
  {

      /* Test if the curves have the same representation and opposite direction.*/

      kstat = 2;
      if (kn1 != kn2 || kk1 != kk2) kstat = 0;
      tdist = s6dist(scoef1,scoef2+(kn2-1)*newdim,newdim);
      if(DNEQUAL(tdist,(double)0.0))
      {
          kstat = 0;
      }
      else
      {
          for (s1=scoef1+newdim,s2=scoef2+(kn2-2)*newdim,ki=1;
		 ki<kn1 && kstat>0; s1+=newdim,s2-=newdim,ki++)
          {
            if (DNEQUAL(s6dist(s1,s2,newdim),tdist))
	    {
	        kstat = 0;
	        break;
	    }
	  }
      }
  }


  /* Difference function made.  */

  *jstat = kstat;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
  s6err("s1956",*jstat,kpos);
  goto out;

  /* Error in in lower level routine. */

 err102: *jstat = kstat;
  s6err("s1956",*jstat,kpos);
  goto out;

  /* Error in input. Conflicting dimensions.  */

  err106 : *jstat = -106;
  s6err("s1956",*jstat,kpos);
  goto out;

 out:

  /* Free space occupied by local array.  */

  if (ssurf != SISL_NULL) freearray(ssurf);
  if (ratpc) freeCurve(ratpc);

  return;
}
