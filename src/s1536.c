/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

#define S1536

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1536(double points[],int im1,int im2,int idim,int ipar,
	   int con1,int con2,int con3,int con4, int order1,
	   int order2,int iopen1,int iopen2,SISLSurf **rsurf,int *jstat)
#else
void s1536(points,im1,im2,idim,ipar,con1,con2,con3,
	   con4,order1, order2,iopen1,iopen2,rsurf,jstat)
     double points[];
     int im1;
     int im2;
     int idim;
     int ipar;
     int con1;
     int con2;
     int con3;
     int con4;
     int order1;
     int order2;
     int iopen1;
     int iopen2;
     SISLSurf **rsurf;
     int *jstat;
#endif
/*
************************************************************************
*
* PURPOSE: To compute a B-spline tensor surface interpolating a set
*          of points.
*
* INPUT:
*          points - Array of dimension idim*im1*im2 containing
*                   the positions of the nodes (using the same ordering
*                   as ecoef in the SISLSurf structure).
*
*          im1    - The number of interpolation points in the
*                   first parameter direction.
*
*          im2    - The number of interpolation points in the
*                   second parameter direction.
*
*          idim   - Dimension of the space we are working in.
*
*          ipar   - Flag showing the desired parametrization to be used:
*                   = 1: Mean accumulated cord-length parameterization.
*                   = 2: Uniform parametrization.
*
*
*                          ^ Second par. direction 
*                          |     
*                          |    (2.)
*                          |-----------|
*                          |           |
*                     (3.) |           | (4.) 
*                          |           |
*                          |           |
*                          |-----------|-> First par. direction
*                               (1.)
*
*          con1      - Additional condition along edge 1:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con2      - Additional condition along edge 2:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con3      - Additional condition along edge 3:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          con4      - Additional condition along edge 4:
*                           = 0: No additional condition.
*                           = 1: Zero curvature.
*
*          order1    - Order of surface in first parameter direction.
*
*          order2    - Order of surface in second parameter direction.
*
*          iopen1    - Open/closed parameter in first parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed surface.
*                      = -1 : closed, periodic surface.
*
*          iopen2    - Open/closed parameter in second parameter direction.
*                      =  1 : open surface.
*                      =  0 : closed surface.
*                      = -1 : closed, periodic surface.
*
*
* Output:
*          rsurf - Pointer to the surf produced
*          jstat  - Status variable
*                    < 0 - Error.
*
* Method:
*     First, a suitable parametrization is calculated according
*     to indicator variable ipar.
*     Then the interpolation is accomplished by using a one dimensional
*     routine for spline interpolation called several times. 
*     First, the datapoints
*     are considered to be idim*im1 dimentional and so on...
*
*
* REFERENCES :
*
* CALLS      : s1536
*
* WRITTEN BY : Christophe Rene Birkeland, SINTEF, May 1993.
*
*********************************************************************
*/                                                               
{
  int kstat=0;        /* Status variable                             */
  int kpos=0;         /* Position of error                           */
  double *par1=SISL_NULL;    /* Transposed positions (in rpos)              */
  double *par2=SISL_NULL;    /* Transposed derivatives (in rder)            */
  
  
  /* Check input */        
  
  if (ipar < 1 || ipar > 3) goto err102;
  
  /* Generate parametrizations */
  
  s1528(idim, im1, im2, points, ipar, iopen1, iopen2, &par1, &par2, &kstat);
  if(kstat < 0) goto error;

  /* Interpolation */

  s1537(points,im1,im2,idim,par1,par2,con1,con2,con3,
	con4,order1, order2,iopen1,iopen2,rsurf,&kstat);
  if(kstat < 0) goto error;

  /* Success */
  
  *jstat = 0;
  goto out;  
  
  /* Error in input data. */

  err102: *jstat = -102;
    s6err("s1536",*jstat,kpos);
    goto out;
    
  /* Error in lower level routine. */

  error:  *jstat =kstat;
    s6err("s1536",*jstat,kpos);
    goto out;
  
  out:
    if(par1 != SISL_NULL) freearray(par1);
    if(par2 != SISL_NULL) freearray(par2);
    return;
}
