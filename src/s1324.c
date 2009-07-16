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
 * $Id: s1324.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1324

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1324(double ecentr[],double aradiu,double enorm[],int idim,
	   double carray[],int *jstat)
#else
void s1324(ecentr,aradiu,enorm,idim,carray,jstat)
     double ecentr[];
     double aradiu;
     double enorm[];
     int    idim;
     double carray[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To make two matrix of dimension 4x4
*              describing a 3-D circle as two implicit functions.
*
*
* INPUT      : ecentr - Center of the circle
*              aradiu - Radius of the circle
*              enorm  - Normal vector of circle plane
*              idim   - The dimension of the space the cirle lies in
*
*
*
* OUTPUT     : carray - The description of the circle. Outside
*                       this function the space for this array must be
*                       allocated. The need is 32 double variables.
*                       First the matrix for the sphere is stored,
*                       then the matrix of the plane.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : The circle is described as an intersection between a
*              cylinder and the plane. The matrix describing the
*              cylinder is put first in the output array, the matrix
*              describing the plane follows then.
*              
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway, 29-June-1988
*                                    
*********************************************************************
*/
{
  int kdimp1;         /* Dimension of matrix kdimp1 = idim + 1         */
  int kdimp2;         /* idim + 2                                      */
  int kstop;          /* Stop condition for for loop                   */
  int ki;             /* Running variables in loop                     */
  int kpos=0;         /* Position of error                             */
  int kstat;          /* Status variable                               */
  
  
  
  /* Test i legal input */
  if (idim != 3) goto err104;
  
  kdimp1 = idim + 1;
  kdimp2 = idim + 2;
  kstop  = 2*kdimp1*kdimp1;
  
  for (ki=0;ki<kstop;ki++)
    {
      carray[ki] = (double)0.0;
    }
  
  /* Make description of cylinder */
  
  s1322(ecentr,enorm,aradiu,idim,1,carray,&kstat);
  if (kstat<0) goto error;
  
  
  /* Make description of plane, element (1,4), (2,4) and (3,4) */
  
  carray[28] = enorm[0];
  carray[29] = enorm[1];
  carray[30] = enorm[2];
  
  /* Make element (4,4) */
  
  carray[31] = -s6scpr(enorm,ecentr,idim);
  
  *jstat = 0;
  goto out;
  
  /* Dimension not 3 */
 err104: *jstat = -104;
  s6err("s1324",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine */
 error: *jstat = kstat;
  goto out;
  
  
 out:
  return;
}
