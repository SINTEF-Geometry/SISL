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
 * $Id: s1435.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S1435

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1435(SISLSurf *ps1,int iedge,SISLCurve **rcedge,double *cpar,int *jstat)
#else
void s1435(ps1,iedge,rcedge,cpar,jstat)
     SISLSurf   *ps1;
     int    iedge;
     SISLCurve  **rcedge;
     double *cpar;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Pick given edge-curve of B-spline surface.
*
*
*
* INPUT      : ps1    - Pointer to surface.
*              iedge  - Number of surface. See figure below.
*
*                           -----------------
*                           !   iedge=2     !
*                           !               !
*                    iedge=3!               !iedge=1
*                           !               !
*                           !               !
*                           !               !
*                           -----------------
*                               iedge=0
*
*
*
* OUTPUT     : rcedge - SISLEdge curve.
*              cpar   - Parameter value of edge in constant direction.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      : s1436 - Pick curve with constant second parameter.
*              s1437 - Pick curve with constant first parameter.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
*
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status parameter.                        */
  int kpos = 0;         /* Position of error.                             */
  double tstart1,tend1; /* Endpoints of parameter interval in first 
			   direction.                                     */
  double tstart2,tend2; /* Endpoints of parameter interval in second 
			   direction.                                     */
  double tpar;          /* Parameter value of curve in constant parameter
			   direction.                                     */
  
  /* Fetch endpoints of parameter intervals.  */
  
  tstart1 = *(ps1->et1 + ps1->ik1 - 1);
  tend1 = *(ps1->et1 + ps1->in1);
  tstart2 = *(ps1->et2 + ps1->ik2 - 1);
  tend2 = *(ps1->et2 + ps1->in2);
  
  /* Find constant parameter of edge. */
  
  if (iedge == 0) tpar = tstart2;
  else if (iedge == 1) tpar = tend1;
  else if (iedge == 2) tpar = tend2;
  else if (iedge == 3) tpar = tstart1;
  
  if (iedge == 0 || iedge == 2)
    {
      
      /* Pick curve with constant second parameter.  */
      
      s1436(ps1,tpar,rcedge,&kstat);
      if (kstat < 0) goto error;
    }
  else if (iedge == 1 || iedge == 3)
    {
      
      /* Pick curve with constant first parameter.  */
      
      s1437(ps1,tpar,rcedge,&kstat);
      if (kstat < 0) goto error;
    }
  
  /* SISLCurve picked.  */
  
  *cpar = tpar;
  *jstat = 0;
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1435",*jstat,kpos);
  goto out;
  
 out: return;
}
