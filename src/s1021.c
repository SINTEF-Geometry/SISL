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
 * $Id: s1021.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1021

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1021(double bottom_pos[], double bottom_axis[], double ellipse_ratio, 
	 double axis_dir[], double height, SISLSurf **cyl, int *stat)
#else
void s1021(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	      height, cyl, stat) 
     double bottom_pos[];
     double bottom_axis[];
     double ellipse_ratio;
     double axis_dir[];
     double height;
     SISLSurf  **cyl;
     int    *stat;      
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To describe a truncated cylinder as a NURBS. The cylinder 
*              can be elliptic.
*             
*
* INPUT      : bottom_pos    - Center point of the bottom
*              bottom_axis   - One of the bottom axis (major or minor)
*              ellipse_ratio - Ratio betwwen the other axis and bottom_axis
*              axis_dir      - Direction of the cylinder axis
*              height        - Height of the cone, can be negative.
*
*
* OUTPUT     : 
*              stat          - status messages  
*                                             > 0      : warning
*                                             = 0      : ok
*                                             < 0      : error
*              cyl           - Pointer to the cylinder produced
*
* METHOD     : The cylinder is made as a rules surface between two NURBS ellipses.
*
*
* REFERENCES :
*
*-                                      
* CALLS      : 
*
* WRITTEN BY : Johannes Kaasa, SI, Oslo, Norway, Jan. 93
*
*********************************************************************
*/
{
   int kstat;                      /* Status variable.   */
   int kpos=0;                     /* Position of error. */
   double cone_angle = (double)0.; /* No conical slope.  */

   s1022(bottom_pos, bottom_axis, ellipse_ratio, axis_dir, 
	     cone_angle, height, cyl, &kstat);
   if (kstat < 0) goto error;
  
   *stat = 0;
   goto out;
  
   /* Error in lower level routine. */
      
   error:
      *stat = kstat;
      s6err("s1021", *stat, kpos);
      goto out;
  
   out:  
      return;
}
    
