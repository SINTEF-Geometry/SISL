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
 * $Id: s1388.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1388

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1388(SISLSurf *ps1,double *gcoons[],int *jnumb1,int *jnumb2,int *jdim,int *jstat)
#else
void s1388(ps1,gcoons,jnumb1,jnumb2,jdim,jstat)
     SISLSurf   *ps1;
     double *gcoons[];
     int    *jnumb1;
     int    *jnumb2;
     int    *jdim;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Convert a B-spline surface of order up to 4 in both
*              direction to a mesh of Coons patches with uniform
*              parametrization. The function assumes that the B-spline
*              surface is C1.
*
*
* INPUT      : ps1    - Pointer to the surface to be converted
*
*
* OUTPUT     : gcoons - Array containing the sequence of Coons patches.
*                       The total number of patches is jnumb1*jnumb2.
*                       The patches are stored in sequence with jdim*16
*                       doubles for each patch. For each corner of the patch
*                       we store in sequence position, derivative in first
*                       direction, derivative in second direction and twist.
*                       This array is allocated inside the routine and must
*                       be released by the calling routine.
*
*              jstat  - status messages  
*                                         = 1      : Orders too high surface
*                                                    interpolated
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : For each polynomial patch the corner position, (1,0)-
*              (0,1)- and (1,1)-derivative is calculated.
*
* USE:         int knumb1,knumb2,kdim,kstat;
*              double *scoons;
*              SISLSurf *qs1;
*               .
*               .
*              s1388(qs1,&scoons,&knumb1,&knumb2,&kdim,&kstat);
*               .
*
*              If one of the orders of the surface is greater than four
*              (*jstat==1) then degree reduction (order reduction) should
*              be applied before using this routine to get a satisfactory
*              representation of the surface by Coons patches.
*              The degree reduction routine is s1348.
*
* REFERENCES :
*
*-
* CALLS      : s1424, s6err
*
* WRITTEN BY : Tor Dokken, SI, Norway, 1988-11
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int kder=1;         /* Calculate all first derivatives                 */
  int klfs=0;         /* Pointer into first knot vector                  */
  int klft=0;         /* Pointer into second knot vector                 */
  int kdumlfs;        /* Temporary pointer into first knot vector        */
  int kdumlft;        /* Temporary pointer into second knot vector       */
  int ksize;          /* Number of doubles to store a Coons patch        */
  
  int kj;          /* Control variables in for loop                   */
  int kn1;            /* Number of vertices in first parameter direction */
  int kn2;            /* Number of vertices in first parameter direction */
  int kk1;            /* Order in first parameter direction              */
  int kk2;            /* Order in first parameter direction              */
  double *st1;        /* Knots in first parameter direction              */
  double *st2;        /* Knots in second parameter direction             */
  double spar[2];     /* Current parameter value                         */
  double sparx[2];    /* Temporary parameter value                       */
  double tdiff1;      /* Length of parameter interval in direction 1     */
  double tdiff2;      /* Length of parameter interval in direction 2     */
  double *scorn1;     /* Pointer to corner 1 of current Coons patch      */
  double *scorn2;     /* Pointer to corner 2 of current Coons patch      */
  double *scorn3;     /* Pointer to corner 3 of current Coons patch      */
  double *scorn4;     /* Pointer to corner 4 of current Coons patch      */
  double tdum;        /* Temporary variable                              */
  
  kn1 = ps1->in1;
  kn2 = ps1->in2;
  kk1 = ps1->ik1;
  kk2 = ps1->ik2;
  kdim = ps1 ->idim;
  st1 = ps1 -> et1;
  st2 = ps1 -> et2;
  
  /* Calculate number of doubles to store a Coons patch */
  
  ksize = kdim*16;
  
  
  /* Allocate array for storage of the coefficients */
  
  *gcoons = newarray((kn1*kn2*16*kdim),DOUBLE);
  if (*gcoons == SISL_NULL) goto err101;
  
  klft = kk2 - 1;
  
  *jnumb2 = 0;
  
  scorn1 = *gcoons;
  
  while (klft < kn2)
    {
      *jnumb1 = 0;
      klfs = kk1 - 1;
      while (klfs < kn1)
        {
	  
	  /* Set pointers to the corners */
	  
	  scorn2 = scorn1 + 4*kdim;
	  scorn3 = scorn2 + 4*kdim;
	  scorn4 = scorn3 + 4*kdim;
	  
	  spar[0] = st1[klfs];
	  spar[1] = st2[klft];
	  
	  /* The parameter pair spar describes the lower left corner of the
	     current polynomial patch. By evaluating at spar we get
	     pointers into the knot vectors telling where we are. These are
	     klfs and klft t. We can calulate in the other corners by
	     just adding one to these pointers and find their parameter value.
	     */
	  
	  /* Calulate first corner of patch */
	  
	  s1424(ps1,kder,kder,spar,&klfs,&klft,scorn1,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Find length of parameter intervals */
	  
	  tdiff1 = st1[klfs+1] - st1[klfs];
	  tdiff2 = st2[klft+1] - st2[klft];
	  
	  /* Calculate second corner of patch */
	  
	  sparx[0] = st1[klfs+1];
	  sparx[1] = spar[1];
	  kdumlfs = klfs;
	  kdumlft = klft;
	  
	  s1424(ps1,kder,kder,sparx,&kdumlfs,&kdumlft,scorn2,&kstat);
	  if (kstat<0) goto error;
	  
	  
	  /* Calculate third corner of patch */
	  
	  sparx[0] = spar[0];
	  sparx[1] = st2[klft+1];
	  kdumlfs = klfs;
	  kdumlft = klft;
	  
	  s1424(ps1,kder,kder,sparx,&kdumlfs,&kdumlft,scorn3,&kstat);
	  if (kstat<0) goto error;
	  
	  /* Calculate third corner of patch */
	  
	  sparx[0] = st1[klfs+1];
	  sparx[1] = st2[klft+1];
	  kdumlfs = klfs;
	  kdumlft = klft;
	  
	  s1424(ps1,kder,kder,sparx,&kdumlfs,&kdumlft,scorn4,&kstat);
	  if (kstat<0) goto error;                 
	  
	  /* Scale derivatives to match uniform parametrization */
	  
	  
	  /* Scale derivatives in first direction */
	  
	  for (kj=kdim;kj<2*kdim;kj++)
            {
	      scorn1[kj] *= tdiff1;
	      scorn2[kj] *= tdiff1;
	      scorn3[kj] *= tdiff1;
	      scorn4[kj] *= tdiff1;
            }
	  
	  /* Scale derivatives in second direction */
	  
	  for (kj=2*kdim;kj<3*kdim;kj++)
            {
	      scorn1[kj] *= tdiff2;
	      scorn2[kj] *= tdiff2;
	      scorn3[kj] *= tdiff2;
	      scorn4[kj] *= tdiff2;
            }
	  
	  /* Scale twist */
	  
	  tdum = tdiff1*tdiff2;
	  
	  for (kj=3*kdim;kj<4*kdim;kj++)
            {
	      scorn1[kj] *= tdum;
	      scorn2[kj] *= tdum;
	      scorn3[kj] *= tdum;
	      scorn4[kj] *= tdum;
            }
	  
	  /* Update position of first point of Coons patch */
	  
	  scorn1 += ksize;
	  
	  klfs += 1;
	  *jnumb1 +=1;
        }
      
      klft += 1;
      *jnumb2 +=1;
    }
  
  /* The array is probably too big for the Coons patches, decrease the
     array size */
  
  
  /* Allocate array for storage of the coefficients */
  
  *gcoons = increasearray(*gcoons,((*jnumb1)*(*jnumb2)*16*kdim),DOUBLE);
  if (*gcoons == SISL_NULL) goto err101;
  
  
  /* Test if order to high */
  
  *jdim = kdim;
  
  if (kk1>4 || kk2 > 4) goto war01;
  
  *jstat = 0;
  
  goto out;
  
  /* Orders too high */
  
 war01:*jstat=1;
  goto out;
  
  
  
  
  /* Error in scratch allocation */
  
 err101: *jstat = -101;
  s6err("s1388",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1388",*jstat,kpos);
  goto freeout;
  
  /* Some error has occured free allocated space */
  
 freeout:
  if (*gcoons != SISL_NULL) freearray(*gcoons);
  goto out;
  
 out:
  return;
}

