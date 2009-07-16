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
 * $Id: s1957.c,v 1.2 2001-03-19 15:58:57 afr Exp $
 *
 */


#define S1957

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1957(SISLCurve *pcurve,double epoint[],int idim,double aepsco,double aepsge,
           double *gpar,double *dist,int *jstat)
#else
void s1957(pcurve,epoint,idim,aepsco,aepsge,gpar,dist,jstat)
     SISLCurve    *pcurve;
     double   epoint[];
     int      idim;
     double   aepsco;
     double   aepsge;
     double   *gpar;
     double   *dist;
     int      *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find the closest point between the curve pcurve and the
*              point epoint. The method is fast and should work well
*              in clear cut cases but does not guarantee finding
*              the right solution. As long as it doesn't fail,
*              it will find exactly one point.
*
*
*
* INPUT      : pcurve - Pointer to the curve in the closest point problem.
*              epoint - The point in the closest point problem.
*              idim   - Dimension of the space in which epoint lies.
*              aepsco - Computational resolution.
*              aepsge - Geometry resolution.
*
*
*
* OUTPUT     : gpar   - The parameter value of the closest point
*                       in the parameter interval 
*                       of the curve.
*              dist   - The closest distance between curve and point.
*              jstat  - status messages  
*                                         > 0      : warning
*                                         = 0      : Point found by iteration
*                                         = 1      : Point lies at an end
*                                         < 0      : error
*
*
* METHOD     : Find an initial guess solution by finding, essentially,
*              the closest control point to the given point
*              and estimating the corresponding parameter value, s1959.
*              This value is then the starting point for a Newton
*              iteration in s1771. The distance of this solution
*              is then compared with the distance of the end points
*              (from the given point) and the minimum is returned.
*              
*              
*
*
* REFERENCES :
*
*- 
* CALLS      : s1959,s1771,newPoint,freePoint.
*
* WRITTEN BY : Michael Floater, SI, 91-10.
*
*********************************************************************
*/                                                               
{                                                                     
  double dist1,dist2;       /* Distances of endpoints from epoint.       */
  double cldist;            /* Current minimum distance.                 */
  double endpt[3];          /* Coeffs of an end point of the curve.      */
  double clspt[3];          /* Coeffs of closest point of the curve.     */
  double enext;             /* Initial guess for iteration               */
  double estart,eend;       /* Parameter area for Newton iteration.      */
  double *et=SISL_NULL;          /* Knot vector.                              */
  int ik;                   /* Order of curve.                           */
  int in;                   /* Number of control points of curve.        */
  int kleft=0;              /* Dummy used in evualation.                 */
  int kstat = 0;            /* Local status variable.                    */
  int kpos = 0;             /* Error position.                           */
  double gpos;              /* Parameter of closest point on curve.      */
  double clgpar;            /* Parameter of current closest point.       */
  SISLPoint *ppoint = SISL_NULL; /* epoint in SISLPoint form.                 */
  
  /* Test input.  */
  
  if (idim != 2 && idim != 3) goto err105;
  if (pcurve->idim != idim) goto err106;

  /* Set up local variables. */

  et = pcurve->et;
  ik = pcurve->ik;
  in = pcurve->in;


  /* Evaluate the curve at its end points. */
  
  s1221(pcurve,0,et[ik-1],&kleft,endpt,&kstat);
  if (kstat < 0) goto error;
  
  dist1 = s6dist(epoint,endpt,idim);
  
  s1221(pcurve,0,et[in],&kleft,endpt,&kstat);
  if (kstat < 0) goto error;
  
  dist2 = s6dist(epoint,endpt,idim);
  
  /* Find the closest end point and store the result. */

  *jstat = 1;

  if(dist1 < dist2)
  {
      cldist = dist1;
      clgpar = et[ik-1];
  }
  else
  {
      cldist = dist2;
      clgpar = et[in];
  }

  
  /* Now try the interior of the curve. */

  ppoint = newPoint(epoint,idim,1);
  if(ppoint == SISL_NULL) goto err101;

  /* Find a good guess point based on finding the closest control
     point and its corresponding parameter values. */

  s1959(ppoint,pcurve,&enext,&kstat);
  if(kstat < 0) goto error;


  /* Do the Newton iteration. */
    
  estart=et[ik-1];
  eend=et[in];

  s1771(ppoint,pcurve,aepsge,estart,eend,enext,&gpos,&kstat);
  if(kstat >= 0)
  {
      /* Closest point found. */
      /* Find distance to compare with end points. */
    
      s1221(pcurve,0,gpos,&kleft,clspt,&kstat);
      if (kstat < 0) goto error;
      
      dist1 = s6dist(epoint,clspt,idim);

      if(dist1 < cldist)
      {
          cldist = dist1;
          clgpar = gpos;
          *jstat = 0;
      }
  }

  /* Return the result. */

  *gpar = clgpar;
  *dist = cldist;


  goto out;
  
  /* Error in space allocation.  */
  
 err101: *jstat = -101;
  s6err("s1957",*jstat,kpos);
  goto out;
  
  /* Error in input. Dimension not equal to 2 or 3.  */
  
 err105: *jstat = -105;
  s6err("s1957",*jstat,kpos);
  goto out;
  
  /* Dimensions conflicting.  */
  
 err106: *jstat = -106;
  s6err("s1957",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine.  */
  
  error : *jstat = kstat;
  s6err("s1957",*jstat,kpos);
  goto out;
  
 out:
  
  /* Free allocated space.  */
  
  if (ppoint != SISL_NULL) freePoint(ppoint);
  
  return;
}                                               
                                           
                       

