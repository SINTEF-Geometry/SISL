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
 * $Id: s1792.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1792

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double 
s1792(double et[],int ik,int in)
#else
double s1792(et,ik,in)
     double et[];
     int    ik;
     int    in;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Finding a good subdividing parametric value.
*
*
*
* INPUT      : et[]     - The knot vector.
*              ik       - The order of the spline function.
*              in       - The number of basic spline function.
*
*
*
* OUTPUT     : s1792    - The subdivision value.
*
*
* METHOD     : This function is written to match s1791().
*
*
* REFERENCES :
*
*
* WRITTEN BY : Arne Laksaa, SI, 89-04.
*
*********************************************************************
*/                                     
{
  if (in > ik)
    {
      int kpar = (in + ik)/2;
      
      if (DNEQUAL(et[ik-1],et[kpar]) || DNEQUAL(et[in],et[kpar]))
        return  et[kpar];
    }
  
  return (et[ik-1]+et[in])*(double)0.5;
}
