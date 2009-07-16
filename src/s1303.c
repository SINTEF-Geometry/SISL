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
 * $Id: s1303.c,v 1.4 1994-08-23 09:18:08 pfu Exp $
 *
 */
#define S1303

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1303(double epstrt[],double aepsge,double angle,double epcent[],
	   double eaxis[],int idim,SISLCurve **rc,int *jstat)
#else
void s1303(epstrt,aepsge,angle,epcent,eaxis,idim,rc,jstat)
     double epstrt[];
     double aepsge;
     double angle;
     double epcent[];
     double eaxis[];
     int    idim;
     SISLCurve  **rc;
     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To create a B-spline approximating a circular arc
*              around the axis defined by the center point epcent[],
*              the axis eaxis[] a start point epstrt[] and rotational
*              angle. The maximal deviation between the true
*              circular arc and the approximation to the arc
*              is controlled by aepsge.
*
*
* INPUT      : epstrt - Start point of circular arc
*              aepsge - Maximal deviation allowed between true circular arc
*                       and generated curve.
*              angle  - The rotational angle. Counter clockwise around axis.
*                       If the rotational angle is outside <-PI,+PI> then
*                       a closed curve is produced.
*              epcent - SISLPoint of axis of circle
*              eaxis  - Normal vector to plane of circle
*              idim   - The dimension of the space (=3)
*
* OUTPUT     :
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*              rc     - Pointer to circle approximation produced.
*                       NB! rc is declared as point to a curve pointer!
*                           Use & in front of actual parameter if it
*                           declared as pointer to curve.
*
* METHOD     : First the maximal distance between the curve and the
*              rotational axis is determined. Then by comparing this
*              with aepsge the allowed relative error is found. This
*              relative error and the rotational angle is used for
*              generating a normalized circle segment spanning the
*              actual angle. This circle is then translated to generate
*              the actual rows of control vertices of the surface
*
* EXAMPLE OF USE:
*              SISLCurve *qr;
*              int    kstat;
*              .
*              .
*              s1303(ep,aepsge,angle,epcnet,eaxis,idim,&qr,&kstat);
*
* REFERENCES :
*
*-
* CALLS      : s6norm, s6scpr, s1301,s6rotax, s6mvec, s6err
*
*
* WRITTEN BY : Tor Dokken, SI, Oslo, Norway. 24. May 1988
* REVISED BY : J. Kaasa, SI, Aug. 92 (Made a proper handling of 2D circles).
* Revised by : Paal Fugelli, SINTEF Oslo, Norway, 23/08-1994. Added handling of
*              kstat after call to s1301().
*
*********************************************************************
*/
{
  double *scoef1;         /* Pointer to vertices of circle segment      */
  int    kn1;             /* Number of vertice of circle segment        */
  double sdiff[3];        /* Array for storing differences              */
  double saxis[3];        /* Array for storing normalized eaxis         */
  int    kj;              /* Control variable in loop                   */
  double tlength;         /* Variable used for length calculation       */
  double treler;          /* Variable used for relative error           */
  double smat[16];        /* Transformation matrix                      */
  int    kstat;           /* Status variable                            */
  double tfak;            /* Value of cross product                     */
  int    kpos = 1;        /* Possition of error                         */

  double tdum;            /* Intermediate storing facility              */
  int    klimit;          /* Limit in for loop                          */


  /* The routine is only working for dimension 2 & 3 */

  if (idim != 2 && idim != 3) goto err104;

  if (idim == 2)
  {
     /* Make difference between start point and point on axis */

     s6diff(epstrt, epcent, idim, sdiff);

     /* Find length of this vector, which is the radius of the circle */

     tlength = s6length(sdiff,idim,&kstat);
  }

  else if (idim == 3)
  {
     /* Normalize axis direction */

     (void)s6norm(eaxis,idim,saxis,&kstat);

     /* Find distance between axis and start point of circle */

     /* Make difference between start point and point on axis */

     s6diff(epstrt, epcent, idim, sdiff);
     tfak = s6scpr(sdiff,saxis,idim);

     /* Find vector normal to axis going to vertex by subtracing the
     component of the difference vector along the axis */

     for (kj=0;kj<3;kj++)
       {
         sdiff[kj] = sdiff[kj] - tfak*saxis[kj];
       }

     /* Find length of this vector, which is the radius of the circle */

     tlength = s6length(sdiff,idim,&kstat);
  }

  /* Calculate relative error. If relative error  <= 0, the whole curve
  lies on the axis  */

  if (tlength <= (double)0.0) goto err127;
  treler = aepsge/tlength;


  /* Calculate normalized circle */

  s1301(treler,angle,idim,rc,&kstat);
  if (kstat < 0)  goto error;

  /* Make local variables for curve description */

  scoef1 = (*rc) -> ecoef;
  kn1    = (*rc) -> in;


  if (idim == 2)
  {

     /* Transform the circle approx into the right position. */

     klimit = idim*kn1;
     for (kj = 0; kj < klimit; kj += 2)
     {
	tdum = epcent[0] + scoef1[kj]*sdiff[0] - scoef1[kj + 1]*sdiff[1];
	scoef1[kj + 1] = epcent[1] + scoef1[kj + 1]*sdiff[0]
	   + scoef1[kj]*sdiff[1];
	scoef1[kj] = tdum;
     }
  }

  else if (idim == 3)
  {

     /* Transform the circle approx into the right position, first make
     transformation matrix */

     s6rotax(epcent,saxis,epstrt,smat,&kstat);

     /* Transform the vertices into right position */

     s6mvec(smat,scoef1,kn1,scoef1);
  }

  *jstat = 0;
  goto out;


  /* Error in input, dimension not equal to 3 */

 err104:
  *jstat = -104;
  s6err("s1303",*jstat,kpos);
  goto out;

  /* Error in input, point lies on axis */

 err127:
  *jstat = -127;
  s6err("s1303",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */

 error:
  *jstat = kstat;
  s6err("s1303", *jstat, kpos);
  goto out;

 out:
  return;
}
