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
 * $Id: s1531.c,v 1.3 2005-02-28 09:04:48 afr Exp $
 *
 */


#define S1531

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1531(double ea[],int idim,int in1,int in2,double **eb,int *jstat)
#else
void s1531(ea,idim,in1,in2,eb,jstat)
     double ea[];
     int idim;
     int in1;
     int in2;
     double **eb;
     int    *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute the transpose in the last two indices, of the
*          matrix given by ea and to output it as eb.
*
* INPUT:
*          ea     - Array of dimension idim*in1*in2 containing the
*                   matrix to be transformed.
*          idim   - The length of the first index of ea.
*          in1    - The length of the second index of ea.
*          in2    - The length of the third index of ea.
*
* OUTPUT:
*          eb     - Pointer to the output array containing the
*                   transposed matrix (dimension idim*in1*in2).
*          jstat  - Status variable
*                    < 0 - Memory allocation error.
* METHOD:
*
* REFERENCES :
*
* CALLS:
*
* WRITTEN BY: Michael Floater, SI, June 1992.
*
*********************************************************************
*/
{
  int i,j,jbase;       /* Loop variable                             */
  int ki,kj,kk;        /* Loop variable                             */
  int idiff;           /*                                           */
  int kpos=0;          /* Position of error                         */
  double *mat=SISL_NULL;    /* Temporary output matrix                   */


  mat = newarray(idim*in1*in2, DOUBLE);
  if(mat == SISL_NULL) goto err101;

  i = 0;
  jbase = 0;
  idiff = (in1 - 1) * idim;

  for(ki=0; ki<in1; ki++,jbase+=idim)
  {
      for(kj=0,j=jbase; kj<in2; kj++,j+=idiff)
      {
          for(kk=0; kk<idim; kk++,i++,j++)
          {
	      mat[i] = ea[j];
          }
      }
  }

  (*eb) = mat;

  /* Calculation completed */

  *jstat = 0;
  goto out;


  /* Error in space allocation */
 err101: *jstat = -101;
  s6err("s1531",*jstat,kpos);
  goto out;


 out:

  return;
}
