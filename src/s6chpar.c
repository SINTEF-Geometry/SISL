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
 * $Id: s6chpar.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6CHPAR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6chpar(double ecoef1[],int in1,int in2,int idim,double ecoef2[])
#else
void s6chpar(ecoef1,in1,in2,idim,ecoef2)
     double ecoef1[];
     int    in1;
     int    in2;
     int    idim;
     double ecoef2[];
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Change parameter directions of vertices of surface.
*
*
*
* INPUT      : ecoef1 - Vertices of original surface.
*              in1    - Number of vertices in first parameter direction.
*              in2    - Number of vertices in second parameter direction.
*              idim   - Dimension of the space in which the surfac lies.
*
*
* OUTPUT     : ecoef2 - Vertices after the changing of parameter directions.
*
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-11.
*
*********************************************************************
*/
{
  register int ki,kj,kk;  /* Counters.  */
  
  for (ki=0; ki<in1; ki++)
    for (kj=0; kj<in2; kj++)
      for (kk=0; kk<idim; kk++)
	ecoef2[(ki*in2+kj)*idim+kk] = ecoef1[(kj*in1+ki)*idim+kk];
}
