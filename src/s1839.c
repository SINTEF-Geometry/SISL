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
 * $Id: s1839.c,v 1.2 2001-03-19 15:58:54 afr Exp $
 *
 */


#define S1839

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1839(SISLSurf *ps1,double epol[],int in,int idim,int *jstat)
#else
void s1839(ps1,epol,in,idim,jstat)
     SISLSurf   *ps1;
     double epol[];
     int    in;
     int    idim;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform improved box-test in intersection involving
*              a surface. Test the surface against a control-polygon.
*
*
*
* INPUT      : ps1     - Pointer to the surface in the intersection.
*              epol    - The control-polygon in the test.
*              in      - Number of vertices in control-polygon.
*              idim    - The dimension of the space in which the polygon
*                        lies.
*
*
*
* OUTPUT     : *jstat  - status messages  
*                            = 2      : Only edge-intersections possible.
*                            = 1      : Internal intersections possible.
*                            = 0      : No intersections possible.
*                            < 0      : error
*
*
* METHOD     : Pick the diagonal of the surface and the tangents in
*              the corners of the surface if these might be different 
*              from the diagonals. Rotate coordinate system according
*              to the vectors picked and perform a box-test.
*
*
* REFERENCES :
*
*-
* CALLS      : s1834 - Perform rotated box-test.
*
* WRITTEN BY : Tor Dokken, SI, OSlo, Norway.
*
*********************************************************************
*/
{
  int kstat = 0;           /* Local status variable.                        */
  int kpos = 0;            /* Position of error.                            */
  int ki;                  /* Counter.                                      */
  int kn1,kn2;             /* Number of vertices in each parameter
			      direction of surface.                         */
  int kk1,kk2;             /* Order in each parameter direction of surface. */
  int kvec;                /* Number of direction vectors to calculate.     */
  int klap;                /* Indcates whether SISLbox of surface overlap point.*/
  double *scoef;           /* Vertices of surface.                          */
  register double *s1,*s2,
  *s3,*s4;                 /* Pointers used to traverse arrays.             */
  double *sdir = SISL_NULL;     /* Array containing direction vectors.           */
  
  /* Check input.  */
  
  if (idim != 2 && idim != 3) goto err105;
  if (idim != ps1->idim) goto err106;
  
  /* Copy surface to local parameters.  */
  
  kn1 = ps1 -> in1;
  kn2 = ps1 -> in2;
  kk1 = ps1 -> ik1;
  kk2 = ps1 -> ik2;
  scoef = ps1 -> ecoef;
  
  /* Find number of rotations to make.  */
  
  if (kk1 > 2 || kk2 > 2) kvec = 10; else kvec = 2;
  
  /* Allocate space for vectors with which the x-axis is to be parallell. */
  
  sdir = newarray(kvec*idim,double);
  if (sdir == SISL_NULL) goto err101;
  
  /* Make diagonal from lower left to upper right corner of patch.  
     s1 points to the array which contains the results, s3 points to the
     lower left corner and s4 to the upper right corner.                 */
  
  for (s1=sdir,s2=s1+idim,s3=scoef,s4=scoef+idim*(kn1*kn2-1); s1<s2;
       s1++,s3++,s4++)
    *s1 = *s4 - *s3;
  
  /* Make diagonal from upper left to lower right corner of patch. s1
     points to the array which contains the results, s3 points to the
     upper left and s4 to the lower right corner.                       */
  
  for (s1=sdir+idim,s2=s1+idim,s3=scoef+idim*kn1*(kn2-1),
       s4=scoef+idim*(kn1-1); s1<s2; s1++,s3++,s4++)
    *s1 = *s4 - *s3;
  
  if (kvec > 2)
    {
      
      /* The surface is not linear in both parameter directions. Make
	 horizontal and vertical tangent in lower left corner. s1 points
	 to the array which contain the results and s3 to the corner.    */
      
      for (s1=sdir+2*idim,s2=s1+idim,s3=scoef; s1<s2; s1++,s3++)
	{       
	  *s1 = *(s3+idim) - *s3;
	  *(s1+idim) = *(s3+idim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in lower right corner. s1
	 points to the array which contain the results and s3 to the corner.*/
      
      for (s1=sdir+4*idim,s2=s1+idim,s3=scoef+idim*(kn1-1); s1<s2; s1++,s3++)
	{
	  *s1 = *(s3-idim) - *s3;
	  *(s1+idim) = *(s3+idim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in upper left corner. s1
	 points to the result array and s3 to the corner.                  */
      
      for (s1=sdir+6*idim,s2=s1+idim,s3=scoef+idim*kn1*(kn2-1);s1<s2;s1++,s3++)
	{
	  *s1 = *(s3+idim) - *s3;
	  *(s1+idim) = *(s3-idim*kn1) - *s3;
	}
      
      /* Make the horizontal and vertical tangent in upper right corner. 
	 s1 points to the result array and s3 to the corner.             */
      
      for (s1=sdir+8*idim,s2=s1+idim,s3=scoef+idim*(kn1*kn2-1);s1<s2;s1++,s3++)
	{
	  *s1 = *(s3-idim) - *s3;
	  *(s1+idim) = *(s3-kn1*idim) - *s3;
	}
    }
  
  /* Rotate coordinate system according to the vectors found and perform
     box-test. First use the diagonal vectors.                           */
  
  klap = 1;
  s1834(scoef,kn1*kn2,epol,in,idim,sdir,sdir+idim,&kstat);
  if (kstat < 0) goto error;
  klap = kstat;
  
  if (klap == 1)
    {
      s1834(scoef,kn1*kn2,epol,in,idim,sdir+idim,sdir,&kstat);
      if (kstat < 0) goto error;
      klap = kstat;
    }
  
  /* If the box-tests performed till now show overlap and the surface
     is non-linear in at least one direction rotate the geometry according
     to the tangent information gathered.                                 */
  
  ki = 2;
  while (ki<kvec && klap == 1)
    {
      s1834(scoef,kn1*kn2,epol,in,idim,sdir+ki*idim,sdir+(ki+1)*idim,&kstat);
      if (kstat < 0) goto error;
      klap = kstat;
      ki += 2;
    }
  
  /* Improved boxtest performed.  */
  
  *jstat = klap;
  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1839",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2 or 3.  */
  
 err105: *jstat = -105;
  s6err("s1839",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("s1839",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1839",*jstat,kpos);
  goto out;
  
 out: 
  
  /* Free allocated space.  */
  
  if (sdir != SISL_NULL) freearray(sdir);
  
  return;
}
