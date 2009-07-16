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
 * $Id: s6schoen.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */



#define S6SCHOEN

#include "sislP.h"
#if defined(SISLNEEDPROTOTYPES)
double
s6schoen(double et[], int ik, int index)
#else
double s6schoen(et,ik,index)
     double et[];
     int    ik;
     int    index;
#endif

/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To determine the knot value of a specified vertice.
*
*
* INPUT      : et     - Double array of dimension [in+ik] containing
*                       the knot vector.
*              ik     - The polynomial order of the B-splines associated
*                       with et.
*              index  - The vertice index at where the knot values are to 
*                       be computed.
*
*                
*
* INPUT/OUTPUT :
*              s6schoen - The knot value at index.
*
*
* METHOD     : The aim is to calculate the knot value of a vertice, using the
*              Schoenberg spline expression:
*
*               *
*              t  = (t    + .............. +t     )/k-1.
*               i     i+1                    i+k-1
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : Per Evensen,SI, August 1991.
*
*********************************************************************
*/                                     
{
  int i;             /* Loop variable                                   */
  double kval=DZERO; /* knot value variable                             */
  
  for (i=index+1;i<index+ik;i++) kval+=et[i];
  kval = kval/(ik-1);

return(kval);
}

