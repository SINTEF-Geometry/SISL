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
 * $Id: s6dist.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6DIST

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6dist(double epoint1[],double epoint2[],int idim)
#else
double s6dist(epoint1,epoint2,idim)
     double epoint1[];
     double epoint2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the distance between the points epoint1 and
*              epoint2.
*
*
*
* INPUT      : epoint1 - First point in distance calculation.
*              epoint2 - Second point in distance calculation.
*              idim    - Dimension of the space in which the points lie.
*
*
*
* OUTPUT     : s6dist  - Distance between the points.
*
*
* METHOD     : Compute lenght of the vector epoint1-epoint2.
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/                                     
{
  register double *s1,*s2,*s3; /* Pointers used to travers epoint1 and epoint2
				  arrays.                                      */
  register double tdist=DZERO; /* Distance between the points.                 */
  
  for (s1=epoint1,s2=epoint2,s3=epoint1+idim; s1<s3; s1++,s2++)
    tdist += (*s1 - *s2)*(*s1 - *s2);
  
  return(sqrt(tdist));
}
