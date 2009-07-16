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
 * $Id: s6equal.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6EQUAL

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
s6equal(double a1,double a2,double aref)
#else
int s6equal(a1,a2,aref)
     double a1;
     double a2;
     double aref;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Test if two numbers are equal.
*
*
*
* INPUT      : a1     - First number.
*              a2     - Second number.
*              aref   - Reference value.
*
*
*
* OUTPUT     : s6equal - Tells if numbers are equal.
*                        = 1 : Equal.
*                        = 0 : Not equal.
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
  double tval;   /* Number used to test equality.  */
  
  tval = a1 - a2;
  tval += aref;
  tval -= aref;
  
  return(DEQUAL(tval,DZERO));
}


