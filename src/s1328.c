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
 * $Id: s1328.c,v 1.3 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1328

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1328(SISLSurf *psold,double epoint[],double enorm1[],double enorm2[],
	   int idim,SISLSurf **rsnew,int *jstat)
#else
void s1328(psold,epoint,enorm1,enorm2,idim,rsnew,jstat)
     SISLSurf   *psold;
     double epoint[];
     double enorm1[];
     double enorm2[];
     int    idim;
     SISLSurf   **rsnew;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Put the equation of the surface pointed at by psold
*              into two planes given by the point epoint and the normals
*              enorm1 and enorm2.. The result is an equation where the 
*              new two-dimensional surface rsnew is to be equal to origo.
*
*
*
* INPUT      : psold  - Pointer to input surface.
*              epoint - SISLPoint in the planes.
*              enorm1 - Normal to the first plane.
*              enorm2 - Normal to the second plane.
*              idim   - Dimension of the space in which the planes lie.
*
*
*
* OUTPUT     : rsnew  - The new two-dimensional surface.
*              jstat  - status messages  
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
* CALLS      : newSurf   - Create and initialize new surface.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
* REVISED BY : Mike Floater, SI, 91-04.
* DEBUGGED BY : Mike Floater, SI, 94-06. Use scSave.
*
*********************************************************************
*/
{
  int kpos = 0;    /* Position of error.                            */
  int kdim;        /* Dimension of the space in which the output 
		      surface lies.                                 */
  int kn1,kn2;     /* Number of coefficients of surface.            */
  int kk1,kk2;     /* Order of surface.                             */
  int ikind;       /* kind of surface psold is.                     */
  double *scoef = SISL_NULL; /* Coeffecient array of new surface.        */
  double *s1,*s2;  /* Pointers used to traverse scoef.              */
  double *sc=SISL_NULL; /* Pointer used to traverse psold->ecoef.        */
  double *scSave=SISL_NULL; /* Pointer to new vertices in rational case. */
  double *s3;      /* Stop pointer of vertex in psold->ecoef.       */
  double *spoint;  /* Pointer used to traverse the point epoint.    */
  double *snorm1;  /* Pointer used to traverse the normal enorm1.   */
  double *snorm2;  /* Pointer used to traverse the normal enorm2.   */
  double *rscoef;  /* Scaled coefficients if psold is rational      */
  double wmin,wmax;/* min and max values of the weights if rational */
  double scale;    /* factor for scaling weights if rational        */
  int i;           /* loop variable                                 */
  int idimp1;      /* idim+1                                        */
  
  /* Test input.  */
  
  if (idim != psold -> idim) goto err106;
  
  /* Set simple variables of the new surface.  */
  
  kdim = 2;
  kn1 = psold -> in1;
  kn2 = psold -> in2;
  kk1 = psold -> ik1;
  kk2 = psold -> ik2;
  ikind = psold -> ikind;
  
  /* rational surfaces are a special case */
  if(ikind == 2 || ikind == 4)
  {
      /* scale the coeffs so that min. weight * max. weight = 1  */
      idimp1=idim+1;
      rscoef = psold -> rcoef;
      wmin=rscoef[idim];
      wmax=rscoef[idim];
      for(i=idim; i< kn1*kn2*idimp1; i+=idimp1)
      {
          if(rscoef[i] < wmin) wmin=rscoef[i];
          if(rscoef[i] > wmax) wmax=rscoef[i];
      } 
      scale=1.0/sqrt(wmin*wmax);
      if ((sc=newarray(kn1*kn2*idimp1,DOUBLE)) == SISL_NULL) goto err101;
      
      for(i=0; i< kn1*kn2*idimp1; i++)
      {
          sc[i]=rscoef[i]*scale;
      } 

      scSave = sc;
  }
  else
  {
      sc = psold -> ecoef;
  }

  /* Allocate space for coeffecient of the new surface.  */
  
  if ((scoef = newarray(kdim*kn1*kn2,double)) == SISL_NULL) goto err101;
  
  /* Compute coefficients of new surface.  */
  
  for (s1=scoef,s2=s1+kdim*kn1*kn2; s1<s2; s1+=2)
    {
      *s1 = *(s1+1) = 0;
      spoint = epoint;
      snorm1 = enorm1;
      snorm2 = enorm2;
      if(ikind == 2 || ikind == 4)
      {
      /* surface is rational so we're using idim+1 - d homogeneous coords */
          for (s3=sc+idim; sc<s3; sc++,spoint++,snorm1++,snorm2++)
	    {
	      *s1 += ((*s3)*(*spoint) - *sc)*(*snorm1);
	      *(s1+1) += ((*s3)*(*spoint) - *sc)*(*snorm2);
	    }
          sc++;
      }
      else
      {
      /* surface is not rational so we're using ordinary idim - d coords */
          for (s3=sc+idim; sc<s3; sc++,spoint++,snorm1++,snorm2++)
	    {
	      *s1 += (*spoint - *sc)*(*snorm1);
	      *(s1+1) += (*spoint - *sc)*(*snorm2);
	    }
      }
    }
  
  
  if(ikind == 2 || ikind == 4) freearray(scSave);

  /* Create output surface.  */
  
  *rsnew = newSurf(kn1,kn2,kk1,kk2,psold->et1,psold->et2,scoef,1,kdim,1);
  if (*rsnew == SISL_NULL) goto err101;
  
  /* Task done.  */
  
  *jstat = 0;
  goto out;
  
  /* Error in space allocation.  */
  
  err101: *jstat = -101;
    s6err("s1328",*jstat,kpos);
    goto out;
  
  /* Error in input. Confliction dimensions.  */
  
  err106 : *jstat = -106;
    s6err("s1328",*jstat,kpos);
    goto out;
  
  out:
  
  /* Free space allocated for local array.  */
  
    if (scoef != SISL_NULL) freearray(scoef);
    return;    
}
