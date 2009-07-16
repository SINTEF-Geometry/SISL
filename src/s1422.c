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
 * $Id: s1422.c,v 1.2 2001-03-19 15:58:49 afr Exp $
 *
 */


#define S1422

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s1422(SISLSurf *ps1,int ider,int iside1,int iside2,double epar[],int *ilfs,
       int *ilft,double eder[],double enorm[],int *jstat)
#else
void s1422(ps1,ider,iside1,iside2,epar,ilfs,ilft,eder,enorm,jstat)
     SISLSurf   *ps1;
     int    ider;
     int    iside1;
     int    iside2;
     double epar[];
     int    *ilfs;
     int    *ilft;
     double eder[];
     double enorm[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Evaluate the surface pointed at by ps1 at the parameter
*              value epar. Compute ider derivatives from the right hand
*	       or the left hand sides.
*
*
*
* INPUT      : ps1    - Pointer to the surface to evaluate.
*              ider   - Number of derivatives to calculate.
*                       < 0 : No derivative calculated.    
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*	       iside1 - Indicator telling if the derivatives in the first
*			parameter direction is to be calculated from the
*			left or from the right:
*			 <  0 calculate derivative from the left hand side
*			 >= 0 calculate derivative from the right hand side.
*	       iside2 - Indicator telling if the derivatives in the second
*			parameter direction is to be calculated from the
*			left or from the right:
*			 <  0 calculate derivative from the left hand side
*			 >= 0 calculate derivative from the right hand side.
*              epar   - Parameter-value at which to calculate. Dimension
*                       of epar is 2.
*
*                
*
* INPUT/OUTPUT : ilfs  - Pointer to the interval in the knotvector
*                        in first parameter direction where epar[0] 
*                        is found. The relation
*                          
*                          et1[ilfs] <= epar[0] < et1[ilfs+1]
* 
*                        where et1 is the knotvektor should hold.
*                        ilfs is set equal to zero at the first call
*                        to the routine. 
*                ilft  - Corresponding to ilfs in the second parameter
*                        direction.
*             
*
*
*
* OUTPUT     : eder   - Array where the derivative of the curve in
*                       apar is placed. The sequence is position,
*                       first derivative in first parameter direction,
*                       first derivative in second parameter direction,
*                       (2,0) derivative, (1,1) derivative, (0,2) 
*                       derivative, etc. Dimension of eder is 
*                       idim*(1+2+...+(ider+1)).
*              enorm  - Normal of surface. Is calculated if ider >= 1.
*                       Dimension is idim. The normal is not normalized.
*              jstat  - status messages  
*                                         = 2      : Surface is degenerate
*                                                    at the point, normal
*                                                    has zero length
*                                         = 1      : Surface is close to
*                                                    degenerate at the point
*                                                    Angle between tangents,
*                                                    less than angular tolerance
*                                         = 0      : ok
*                                         < 0      : error
*                      
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : 
*
* WRITTEN BY : 
*
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kpos=0;         /* Position of error.                              */
  int kdim;           /* Dimension of the space in which the surface lies. */
  int keder;          /* Integer used in address calculations on eder    */
  int ksp;            /* Integer used in address calculations on sp      */
  int kincre;         /* Increment for address calculations              */
  int ki,kl;          /* Control variables in for loop                   */
  int knumb;          /* Number of elements used for storage of deriv.s  */
  double *sp;         /* Pointer to temporary array                      */
  double sdum[48];    /* Array used in stead of allocation               */
  
  
  /* Allocate array for storage of ider*ider derivatives */
  
  sp = SISL_NULL;
  kdim = ps1 -> idim;
  knumb = kdim*(ider+1)*(ider+1);
  
  /* Only allocate space if sdum is too smaall */
  
  if (knumb>48)
    sp = newarray(knumb,DOUBLE);
  else
    sp = &sdum[0];
  
  if (sp == SISL_NULL) goto err101;
  
  
  /* Evaluate s1422surface.  */
  
  s1425(ps1,ider,ider,iside1,iside2,epar,ilfs,ilft,sp,&kstat);
  
  if (kstat < 0) goto error;
  
  /* Copy required derivatives into eder */
  
  kincre = kdim*ider;
  
  /*  Copy all derivatives of order 0, then of order 1, up to order ider */
  
  for (kl=0,keder=0;kl<=ider;kl++)
    {
      for (ki=0,ksp=kl*kdim ; ki<=kl ; ki++,ksp+=kincre,keder+=kdim)
        {
	  memcopy(eder+keder,sp+ksp,kdim,DOUBLE);
        }
    }
  
  /* Make cross products of tangents, if idim==3 and derivative >0 */
  
  if (ider>0 && kdim ==3)
    {
      double tlen1,tlen2,tnorm,tang=(double)0.0;
      
      s6crss(eder+kdim,eder+2*kdim,enorm);
      
      /*  Make length of tangents and normal */
      
      tlen1 = s6length(eder+kdim,kdim,&kstat);
      tlen2 = s6length(eder+2*kdim,kdim,&kstat);
      tnorm = s6length(enorm,kdim,&kstat);
      
      /*  Calculate angle between tangents */
      
      if (tlen1 != DZERO && tlen2 != DZERO && tnorm != DZERO)
        tang = tnorm/(tlen1*tlen2);
      
      if (tang == DZERO) *jstat = 2;
      else if (tang <= ANGULAR_TOLERANCE) *jstat = 1;   
      else *jstat = 0;
      goto out;
      
    }
  
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1422",*jstat,kpos);
  goto out;
  
 err101: *jstat = -101;
  s6err("s1422",*jstat,kpos);
  
  
 out:
  
  /* Free allocated space (Space only allocated if sdum is too small) */
  
  if (knumb>48)
    if (sp != SISL_NULL) freearray(sp);
  
  return;
}


