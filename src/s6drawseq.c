/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s6drawseq.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6DRAWSEQ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
extern void line3(double,double,double,int);
#else
extern void line3();
#endif

#if defined(SISLNEEDPROTOTYPES)
void 
  s6drawseq(double epoint[],int ipoint)
#else
void s6drawseq(epoint,ipoint)
     double epoint[];
     int ipoint;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Draw a broken line as a sequence of straight lines 
*              described by the array epoint.
*
*
*
* INPUT      : epoint - Array describing the corners of the broken line
*                       to draw. Each corner is given as a 3-dimensional
*                       point. The corners is placed continuous in the
*                       array.
*              ipoint - Number of points in epoint.
*
*
* OUTPUT     : 
*
*
* METHOD     :
*
*
* REFERENCES :
*
*
* REMARK     : This routine is machine-dependant and the internal of 
*              this routine has to be reprogrammed if the function
*              line3 does not exist on the current computer.
*
*
*-
* CALLS      : line3  - Place pen at given position or draw a line from
*                       the current position to the given one.
*
* WRITTEN BY :
*
*********************************************************************
*/
{
  int ki;          /* Counter.                                    */
  double *spoint;  /* Pointer to corner point in the broken line. */
  
  /* Position pen at start of the broken line.  */
  
  line3(epoint[0],epoint[1],epoint[2],0);
  
  /* Draw sequence of line-segments.  */
  
  for (ki=1,spoint=epoint+3; ki<ipoint; ki++,spoint+=3)
    line3(spoint[0],spoint[1],spoint[2],1);
  
  return;
}
                                    
