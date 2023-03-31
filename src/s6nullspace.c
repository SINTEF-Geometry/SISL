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


#define S6NULLSPACE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s6nullspace(double ea[],int im1,int im2,double aepsge,
		 double **nullspace,int *numbvect,int *jstat)
#else
void s6nullspace(ea,im1,im2,aepsge,nullspace,numbvect,jstat)
     double ea[];
     int    im1;
     int    im2;
     double aepsge;
     double **nullspace;
     int    *numbvect;
     int    *jstat;
#endif
/*
***********************************************************************
*
************************************************************************
*
*   PURPOSE : To find the nullspace of a matrix within the tolerance
*             aepsge
*
*
*   INPUT   : ea     - The matrix
*             im1    - The number of rows
*             im2    - The number of columns
*             aepsge - The size of numbers to be regarded zero
*
*
*
*   OUTPUT  : nullspace    - Array of vectors spanning the approximated
*                            nullspace
*             numbvect     - Number of vectors spanning the null space
*             jstat - Status variable.
*                       < 0 : error
*                       = 0 : ok
*                       > 0 : warning
*
*
*-
*   CALLS      :
*
*   WRITTEN BY : Vibeke Skytt, SI, 86-12.
*   Change to nullspace search by: Tor dokken, SINTEF, 94-06
*
************************************************************************
*/
{
  int *ellimvert=NULL; /* Array telling which column has been ellimnated */
  int *ellimhort=NULL; /* Array telling which row has been elliminated */
  int knoneell;        /* Variable for accumulateing none elliminated variables */

  int kpos = 0;   /* Position of error.                               */
  int ki,kj,kl;   /* Counters.                                        */
  int kpek1,kpek2; /* Pointer pair to maximum element */

  double *sellim=NULL; /* The row used as current row for ellimination*/
  double tmax;    /* Variable used for storing maximale values */
  double tscale;  /* Variable used for storing scaling factors */
  double *qs=NULL;/* Temporary array pointer */
  double *qr=NULL;/* Temporary array pointer */
  double temp;    /* Temporary value */
  double *tempmat=NULL; /* Error storing ellimination matrix row by row */
  double *tempres=NULL;

  /* Allocate space for local array.  */

  if ((ellimvert = newarray(im2,int)) == NULL) goto err101;
  if ((ellimhort = newarray(im1,int)) == NULL) goto err101;
  if ((tempmat = newarray(im2*im2,double)) == NULL) goto err101;

  /* Initite ellimvert and ellimhort to -1 */

  for(ki=0;ki<im2;ki++)
    ellimvert[ki]=-1;

  for(ki=0;ki<im1;ki++)
    ellimhort[ki]=-1;

  /* Find largest element in each row.  */

  for (ki=0; ki < MIN(im1,im2-1); ki++)
    {
      /* Elliminate in column with biggest remaining element */

      tmax = (double)0;
      kpek1=-1;
      kpek2=-1;

      for (kj=0,qs=ea; kj<im1; kj++,qs+=im2)
	if(ellimhort[kj]==-1)
	  {
	    /* Row not used as bases for reduction */
	    for(kl=0; kl<im2; kl++)
	      if(ellimvert[kl]==-1)
		{
		  /* column not used in reduction */
		  temp = fabs(qs[kl]);
		  if(temp>=tmax)
		    {
		      tmax=temp;
                      kpek1=kj;
		      kpek2=kl;
		    }
		}
	  }
      /* fprintf(stdout,"\n tmax %lf",tmax); */
      if(tmax>aepsge)
	{
	  /* If kpek1 and kpek2 greater than -1 then an element remain.
	     It tmax <=aepsge then the remaining matrix should be regarded
	     as zero. Elliminate based on row kpek1 in column kpek2 */

	  /* Find the row to use in the ellimination and scale element such
	     that ellement kpek2 is 1 */
	  sellim = &ea[kpek1*im2];
	  tscale = sellim[kpek2];
	  for(kl=0;kl<im2;kl++)
	    sellim[kl]/=tscale;

	  /* Elliminate from all other rows */

	  for (kj=0,qs=ea; kj<im1; kj++,qs+=im2)
	    if(kj != kpek1)
	      {
		/* We exclude ellimination from row used for ellimination */
		tscale = qs[kpek2];
		for(kl=0; kl<im2; kl++)
		  qs[kl]-=tscale*sellim[kl];
	      }
	  ellimhort[kpek1]=1;
	  ellimvert[kpek2]=kpek1;
	}
    }
/* Build matrix */


for(ki=0,knoneell = 0,qs=tempmat;ki<im2;ki++,qs+=im2)
  {
    /* Check i variable elliminated */

    if(ellimvert[ki]==-1)
      {
	/* Variable not elliminated,make row
           that projects current remianing variable
           onto the right position */
	for(kj=0;kj<im2;kj++)
	  qs[kj]=(double)0;
	qs[knoneell]=-(double)1;
        knoneell+=1;
      }
    else
      {
	/* Variable elliminated, ellimvert[ki] gives actual row */
        qr = ea + im2*ellimvert[ki];
	for(kj=0,kl=0;kj<im2;kj++)
	  if(ellimvert[kj]==-1)
	    {
	      qs[kl] = qr[kj];
	      kl++;
	    }
      }
  }


if ((tempres = newarray(knoneell*im2,double)) == NULL) goto err101;

/* The result is stored row wise, with the first knoenell elements
   relevant to the result, store result column wise */

for(ki=0,qs=tempres,qr=tempmat;ki<knoneell;ki++,qs+=im2,qr+=1)
  for(kj=0;kj<im2;kj++)
    qs[kj] = qr[im2*kj];

*nullspace = tempres;
*numbvect = knoneell;

  *jstat = 0;
  goto out;

/* Error in space allocation.  */

err101: *jstat = -101;
        s6err("s6lufacp",*jstat,kpos);
        goto out;

out:

/* Free space occupied by local array.  */

if (ellimhort != NULL) freearray(ellimhort);
if (ellimvert != NULL) freearray(ellimvert);
if (tempmat != NULL)   freearray(tempmat);
return;
}
