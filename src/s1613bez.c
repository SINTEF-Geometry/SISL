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
 * $Id: s1613bez.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1613BEZ

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1613bez(SISLCurve *pc,int idiv,double aepsge,double **gpar,int *jnpar,
	      int *jstat)
#else
void s1613bez(pc,idiv,aepsge,gpar,jnpar,jstat)
     SISLCurve  *pc;
     int idiv;
     double aepsge;
     double **gpar;
     int    *jnpar;
     int    *jstat;
#endif
/*
*********************************************************************
* 
* PURPOSE    : To compute an aproximation of a Bezier curve with a 
*              sequence of lines inside a tolerance of aepsge.
* 
* 
* INPUT      : pc    - Given spline curve.
*              idiv  - Inital number into which to divide the curve.
*              aepsge    - Tolerance.
*
* 
* OUTPUT     : gpar     - Array containing the parameter values of the points.
*                         Array are allocated inside this function,
*              jnpar    - Number of sampling points.
*              jstat     - status messages 
*                          = 2 : warning. Level of recursion too deep.
*                          = 0 : ok 
*                          < 0 : error 
* 
* 
* METHOD     :
*              
*
* REFERENCES : 
*              
*
* USE        : 
*
*-
* CALLS      : newknots, s1018, s6dline, newCurve, s1613bez,
*              s6takeunion, freeCurve
*              
*
* WRITTEN BY : Vibeke Skytt, SI, 08.92.
*
*********************************************************************
*/
{
   static int klevel = 0;  /* Level of recursion.                    */
   int kstat = 0;        /* Local status varaible.                   */
   int ki;               /* Counter.                                 */
   int kord = pc->ik;    /* Order of curve.                          */
   int kdim = pc->idim;  /* Dimension of geometry space.             */
   int knpar;            /* Number of par. values in refinement.     */
   int kbez = 3;         /* Polynomial Bezier curve.                 */
   int kdiv = 10;        /* Initial division of subcurve.            */
   int kparout;          /* Number of divisions of subcurve.         */
   int knunion;          /* Number of parameter values in union vector. */
   int knparts = idiv;   /* Number of division parameter values.     */
   double tstart = *(pc->et+kord-1);   /* Start par. value of curve. */
   double tend = *(pc->et+pc->in);     /* End par. value of curve.   */
   double tint = (tend - tstart)/(double)(idiv+1);  /* Parameter interval
						       between div. points. */
   double tpar;             /* Current parameter values.                */
   double *sparout = SISL_NULL;  /* Array of parameter values from subcurve. */
   double *spar = SISL_NULL;     /* Array of parameter values.        */
   double *sparunion = SISL_NULL;  /* Union of parameter arrays.      */
   double *spar2 = SISL_NULL;    /* Array of parameter values.        */
   double *sta,*stb,*s1;    /* Pointers into coefficient arrays. */
   SISLCurve *qc = SISL_NULL;    /* Refined curve.                    */
   SISLCurve *qbez = SISL_NULL;  /* Bezier segment of refined curve.  */

   /* Check level of recursion.  */
   
   klevel++;
   if (klevel > 200)
   {
      *jstat = 2;
      goto out;
   }
   
   /* Check if the curve is a Bezier curve.  */
   
   if (pc->ik != pc->in) goto err111;

   /* Set up array of division parameter values. First allocate scratch.  */
      
   if ((spar = newarray(idiv,DOUBLE)) == SISL_NULL) goto err101;
     
   for (ki=0, tpar=tstart+tint; ki<idiv; ki++, tpar+=tint)
      spar[ki] = tpar;
   
   /* Fetch the new knots with which to refine the curve. */
   
   newknots(pc->et,pc->in,pc->ik,spar,idiv,REL_PAR_RES,&spar2,&knpar,&kstat);
   if (kstat < 0) goto error;
		  
   /* Refine the curve.  */
		  
   s1018(pc,spar2,knpar,&qc,&kstat);
   if (kstat < 0) goto error;
		  
   /* Traverse Bezier segments and test if linearization is ok, i.e.
      the distance between inner vertices and the line segment between
      the first and last vertex, is less than aepsge.       */
		  
   for (sta=qc->ecoef, ki=0; ki<=idiv; sta=stb+kdim, ki++)
   {
      for (stb=sta+(kord-1)*kdim, s1=sta+kdim; s1<stb; s1+=kdim)
      {
	 if (s6dline(sta,stb,s1,kdim,&kstat) > aepsge) break;
	 if (kstat < 0) goto error;
      }
      if (kstat < 0) goto error;
		     
      if (s1 < stb)
      {
	 /* Precision in linearization not good enough. Divide the
	    Bezier segment and try again. First express the Bezier
	    segment as a curve.   */
	 
	 if ((qbez = newCurve(kord,kord,qc->et+ki*kord,sta,kbez,kdim,0))
	     == SISL_NULL) goto err101;
	 
	 /* Linearize.  */
	 
	 s1613bez(qbez,kdiv,aepsge,&sparout,&kparout,&kstat);
	 if (kstat < 0) goto error;
	 
	 if (kstat == 2)
	 {
	    /* Return to top level.  */
	    
	    if (spar != SISL_NULL) freearray(spar);
	    if (qbez != SISL_NULL) freeCurve(qbez);
	    
	    *jstat = kstat;
	    goto out;
	 }
	 
	 /* Take union of parameter values.    */
	 
	 s6takeunion(spar,knparts,sparout,kparout,&sparunion,&knunion,&kstat);
	 if (kstat < 0) goto error;
	 
	 freearray(spar);
	 spar = sparunion;
	 sparunion = SISL_NULL;
	 knparts = knunion;
	 
	 if (sparout != SISL_NULL) freearray(sparout);
	 sparout = SISL_NULL;
	 
	 freeCurve(qbez);
	 qbez = SISL_NULL;
      }
   }
   
   /* Return parameter values.  */
   
   *gpar = spar;
   *jnpar = knparts;
   
   *jstat = 0;
   goto out;
   
   /* Error in scratch allocation.  */
   
   err101 :
      *jstat = -101;
   goto out;
   
   /* Error. The curve is not of Bezier type.  */
   
   err111 :
      *jstat = -111;
   goto out;
      
   /* Error in lower order routine.  */
   
   error :
      *jstat = kstat;
   goto out;
   
   out:
      /* Reduce level parameter.  */
      
      klevel--;
   
      /* Free scratch used for local arrays.  */
      
      if (spar2 != SISL_NULL) freearray(spar2);
      if (qc != SISL_NULL) freeCurve(qc);
		      
      return;
}
     
