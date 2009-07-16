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
 * $Id: s1791.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1791

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s1791(double et[],int ik,int in)
#else
int s1791(et,ik,in)
     double et[];
     int    ik;
     int    in;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if it is possible to insert new internal nots
*              any further.
*
*
*
* INPUT      : et     - The knot vector.
*              ik     - The order of the curve/surface.
*              in     - The number of the basic functions.
*
*
*
* OUTPUT     : s1791  - Result of the test.
*                       = 0 : It is not possible to insert new knots
*                             any further.
*                       = 1 : It is possible to insert new knots.
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
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/                                     
{
  register double tstart= et[ik - 1];
  register double tend  = et[in];
  register double tmid  = (tstart+tend)*(double)0.5;
  
  /* Check if it is possible to divide the parameter interval.  */
  
  if (DEQUAL(tmid,tstart) || DEQUAL(tmid,tend)) 
    return  0;
  else 
    return  1;
}





