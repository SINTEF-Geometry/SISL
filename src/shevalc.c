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
 * $Id: shevalc.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SHEVALC

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
shevalc(SISLCurve *pc1,int ider,double ax,double aepsge,int *ileft,
	     double eder[],int *jstat)
#else
void shevalc(pc1,ider,ax,aepsge,ileft,eder,jstat)
     SISLCurve *pc1;
     int ider;
     double ax;
     double aepsge;
     int *ileft;
     double eder[];
     int *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To compute the value and ider first derivatives of the
*              B-spline curve pointed to by pc1, at the point with
*              parameter value ax. Use the filtered coefficients of the
*              curve.
*
*
*
* INPUT      : pc1    - Pointer to the curve for which position
*                       and derivatives are to be computed.
*              ider   - The number of derivatives to compute.
*                       < 0 : Error.
*                       = 0 : Compute position.
*                       = 1 : Compute position and first derivative.
*                       etc.
*              ax     - The parameter value at which to compute
*                       position and derivatives.
*              aepsge - Geometry resolution.
*
*                
*
* INPUT/OUTPUT : ileft - Pointer to the interval in the knot vector
*                        where ax is located. If et is the knot vector,
*                        the relation
*                          
*                          et[ileft] <= ax < et[ileft+1]
* 
*                        should hold. (If ax == et[in] then ileft should
*                        be in-1. Here in is the number of B-spline
*                        coefficients.)
*                        If ileft does not have the right value upon
*                        entry to the routine, its value will be changed
*                        to the value satisfying the above condition.
*
*
*
* OUTPUT     : eder   - Double array of dimension [(ider+1)*idim]
*                       containing the position and derivative vectors.
*                       (idim is the number of components of each B-spline
*                       coefficient, i.e. the dimension of the Euclidean
*                       space in which the curve lies.)
*                       These vectors are stored in the following order:
*                       First the idim components of the position vector,
*                       then the idim components of the tangent vector,
*                       then the idim components of the second derivative
*                       vector, and so on.
*                       (The C declaration of eder as a two dimensional array
*                       would therefore be eder[ider+1,idim].)
*              jstat  - Status messages  
*                                         > 0      : Warning.
*                                         = 0      : Ok.
*                                         < 0      : Error.
*
*
* METHOD     : 
*
* REFERENCES :
*
*-
* CALLS      : s1221    - Evaluate curve.
*              s1991    - Make the direction cone of a curve.
*              newCurve - Create new curve object.
*              freeCurve - Free scratch occupied by curve object. 
*
* WRITTEN BY :  Vibeke Skytt, SI, 04.91.
* CORRECTED BY: UJK, SI, 06.91
*********************************************************************
*/                                     
{
  int kstat=0;        /* Local status variable.                          */
  int kdim = pc1->idim;  /* Dimension of geometry space.                 */
  double *scoef=NULL;    /* Array storing filtered coefficients.         */
  double *s1,*s2,*s3,*s4; /* Pointers into coefficient arrays.           */
  SISLCurve *qc = NULL;   /* Curve to evaluate.                          */
  
  /* Make sure that the filtered coefficients of the curve exist.  */
  
  if (kdim == 1)
  {
     /* Create filtered coefficients. */
     
     if ((scoef = newarray(pc1->in,DOUBLE)) == NULL) goto err101;
     
     for (s1=pc1->ecoef, s2=scoef, s3=s1+pc1->in; s1<s3; s1=s4)
     {
	*s2 = *s1;
	for (s2++, s4=s1+1; s4<s3; s4++, s2++)
	{
	   if (fabs((*s4)-(*s1)) < aepsge) *s2 = *s1;
	   else break;
	}
     }
     
     /* Create curve object.  */
     
     if ((qc = newCurve(pc1->in,pc1->ik,pc1->et,scoef,pc1->ikind,
			kdim,0)) == NULL) goto err101;
  }
  else
  {
     /* Check if the cone of the curve exist. */
     
     if (pc1->pdir == NULL)
     {
	/* Create cone.  */
	
	s1991(pc1,aepsge,&kstat);
	if (kstat < 0) goto error;
     }
     
     /* Create curve object.  */
     
     if ((qc = newCurve(pc1->in,pc1->ik,pc1->et,pc1->pdir->esmooth,
			pc1->ikind,kdim,0)) == NULL) goto err101;
  }
  
  /* Evaluate curve.  */
  
  s1221(qc,ider,ax,ileft,eder,&kstat);
  if (kstat < 0) goto error;
	
  /* UJK Let's have a normal exit possibility !*/
  *jstat = 0;
  goto out;

  /* Error in scratch allocation.  */
  err101: *jstat = -101;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  goto out;

out: 
   /* Free scratch occupied by local objects. */
   
   if (scoef != NULL) freearray(scoef);
   if (qc != NULL) freeCurve(qc);
   
   return;
}
