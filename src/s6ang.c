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
 * $Id: s6ang.c,v 1.2 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S6ANG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s6ang(double evec1[],double evec2[],int idim)
#else
double s6ang(evec1,evec2,idim)
     double evec1[];
     double evec2[];
     int    idim;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Compute the angle (in radians) between two vectors
*
*
*
* INPUT      : evec1   - First vector 
*              evec2   - Second vector 
*              idim    - Dimension of the space in which the vectors lie.
*
*
*
* OUTPUT     : s6ang   - Angle in radians between vectors
*
*
* METHOD     : Make cosine of the angle by computing the scalar product,
*              then divide by the length of the two vectors.
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Tor Dokken SI, 88-07.
*              Arne Laksaa SI, 89-07.
*
*********************************************************************
*/                                     
{
  double tscpr,tang,tlength1,tlength2,tcos;
  int    kstat1,kstat2;
  
  tscpr = s6scpr(evec1,evec2,idim);
  
  tlength1 = s6length(evec1,idim,&kstat1);
  tlength2 = s6length(evec2,idim,&kstat2);
  
  if (!kstat1 || !kstat2)
    tang = DZERO;
  else
    {
      tcos = fabs(tscpr/(tlength1*tlength2));
      tcos = MIN((double)1.0,tcos);
      tang = acos(tcos);
    }
  
  return(tang);
}
