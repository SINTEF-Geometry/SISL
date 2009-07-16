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
 * $Id: s6affdist.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6AFFDIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
     s6affdist(double e1[],double e2[],double emat[],int idim)
#else
double s6affdist(e1,e2,emat,idim)
   double e1[];
   double e2[];
   double emat[];
   int    idim;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between two points using an
*              affine metric described by the matrix emat.
*
*
*
* INPUT      : e1     - First point.
*              e2     - Second point.
*              emat   - Matrix of affine metric.
*              idim   - Dimension of geometry space.
*              
*
* OUTPUT     : s6affdist - Distance between two points.
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Vibeke Skytt, SI, 91-03.
*
*********************************************************************
*/
{
   int ki,kj;              /* Counters.  */
   double tdist = DZERO;   /* Distance.  */
   
   for (ki=0; ki<idim; ki++)
      for (kj=0; kj<idim; kj++)
	 tdist += emat[ki*idim+kj]*(e1[ki]-e2[ki])*(e1[kj]-e2[kj]);
   
   tdist = sqrt(idim*tdist);
   
   return tdist;
}
