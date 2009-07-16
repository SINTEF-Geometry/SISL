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
 * $Id: s9conmarch.c,v 1.2 2001-03-19 15:59:02 afr Exp $
 *
 */


#define S9CONMARCH

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s9conmarch(SISLSurf *ps,double alevel,double epar[],int ndir[],int ipoint,
		double *gpar[],int *mpar[],int *jpoint,int *jstat)
#else
void s9conmarch(ps,alevel,epar,ndir,ipoint,gpar,mpar,jpoint,jstat)
     SISLSurf   *ps;
     double alevel;
     double epar[];
     int    ndir[];
     int    ipoint;
     double *gpar[];
     int    *mpar[];
     int    *jpoint;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To check which of the points in the epar array can be
*              connected by marching. If the connection should be
*              to internal point these points are returned. All points
*              are assumed to lie on the boundary or close to the boundary
*              of the patch (DEQUAL says that they lie on the boundary).
*
*              When ipoint = 2, we connect always if marching succeeds!!!!!!!
*              without testing for equality between marching results and
*              given point no 2.
*
* INPUT      : ps        - The surface in intersection.
*              alevel    - The constant value the surface is intersected with.
*              epar[2,*] - Parameter values for the  intersection points along
*                          the boundary.
*              ndir[*]   - Description of the points in epar.
*                          -1 : tangent pointing out
*                           0 : tangent parallel to boundary
*                           1 : tangent pointing in
*                           2 : singular point
*                          The algorithm assumes that first all points
*                          with status -1 and 1 lies first in epar.
*              ipoint    - Number of parameter values
*
*
* OUTPUT     : gpar      - Pointer to parameter values produced tor singular 
*                          points. Is allocated inside the function, must
*                          be released by the calling function.
*              mpar      - Array describing how the input parameter pairs were
*                          connected to other in the input vector and to the
*                          new parameter pairs made in gpar. Is allocated
*                          inside the function, must be released by the
*                          calling function.
*                           -3             - point on closed curve
*                           -2             - internal point onc other curve
*                           -1             - Can not march from point
*                           POSITIVE VALUE - SISLPoint connected to point in epar
*              jpoint    - Number of new parameter pairs produced
*              jstat   -  status messages  
*                                = 2   : Only singular points, not connected
*                                = 1   : Points connected
*                                = 0   : Points not connected..
*                                < 0   : Error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. August  1989
*
*********************************************************************
*/
{
  int kstat;            /* Status variable                             */
  int kpos=0;           /* Position of error                           */
  int *lpar = SISL_NULL;     /* Pointer to output integer array             */
  int ki,kj;
  int kn1,kn2,kk1,kk2;  /* Surface attributes.           */
  double tstart1,tstart2,tend1,tend2; /* Surface attributes.           */
  int ksucc;            /* Success indicator                           */
  double tepsge=1.0;    /* Not used                                    */
  double *spar=SISL_NULL;    /* Pointer to output real array                */
  double scand1[2];     /* Result of iteration process                 */
  double scand2[2];     /* Result of iteration process                 */
  double *sp,*sq;       /* Pointer used in loop                        */
  double tdum1;         /* Max knot value used in DEQUAL comparing.    */
  double tdum2;         /* Max knot value used in DEQUAL comparing.    */

  /* Init */
  kn1 = ps->in1;
  kn2 = ps->in2;
  kk1 = ps->ik1;
  kk2 = ps->ik2;

  tstart1 = ps->et1[kk1-1];
  tend1   = ps->et1[kn1];
  tstart2 = ps->et2[kk2-1];
  tend2   = ps->et2[kn2];

  tdum1 = (double)2.0*max(fabs(tstart1),fabs(tend1));
  tdum2 = (double)2.0*max(fabs(tstart2),fabs(tend2));

  /* Allocate output arrays */
  
  if ((*mpar=newarray(3*ipoint,INT     )) == SISL_NULL) goto err101;
  if ((*gpar=newarray(6*ipoint,DOUBLE)) == SISL_NULL) goto err101;
  
  lpar = *mpar;
  spar = *gpar;
  
  memcopy(spar,epar,2*ipoint,DOUBLE);
  *jpoint = ipoint;
  
  /* Initiate output integer array to point to no points */
  
  for (ki=0 ; ki< 3*ipoint ; ki++) *(lpar+ki) = 0;
  

  /* Loop for all input points. */      
  for (ki=0, sp=spar ; ki< ipoint-1 ; ki++, sp+=2)
    {
      /* Start marching from point ki */

      /* Exclude points already connected and parallell points. */
      if (lpar[ki] != 0 || ndir[ki] == 0) continue;
	  
      /* SISLPoint not marched to */
	  
      s1787(ps,alevel,tepsge,sp,scand1,scand2,&kstat);
      if (kstat<0) goto error;
      if (kstat==0) goto war00;
	  
      /* Run through remaining points to find if scand2 matches any
	 of them. If we've got only two points, we connect them.*/
	  
      ksucc = 0;
	  
      for (kj=ki+1,sq=spar+2*ki+2 ; kj<ipoint ; kj++,sq+=2)
	{
	      
	  /* SISLPoint found */
	      
	  if (DEQUAL(sq[0]+tdum1,scand2[0]+tdum1) && 
	      DEQUAL(sq[1]+tdum2,scand2[1]+tdum2))
	    {
	      /* Accepted end point found */
	      
	      lpar[ki] = kj+1;
	      lpar[kj] = ki+1;
	      ksucc = 1;
	      break;
	    }
	}
      /* If ksucc==0 then one of the searches was not successful */
	  
      if (ksucc==0) goto war00;   

    }
  
  goto success;

 success: 
  *jstat = 1;
  goto out;

  /* No success */
 war00: 
  *jstat=0;
  /* If we got only singular points, set status. */
  if (ndir[0] == 2) *jstat = 2;

  goto out;

  /* Error in space allocation */
 err101: 
  *jstat = -101;
  s6err("s9conmarch",*jstat,kpos);
  goto out;

  /* Error in lower level function */
 error:
  *jstat = kstat;
  s6err("s9conmarch",*jstat,kpos);
  goto out;

 out:;
}
