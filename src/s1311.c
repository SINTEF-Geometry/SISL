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
 * $Id: s1311.c,v 1.3 2001-03-19 15:58:44 afr Exp $
 *
 */


#define S1311

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
double
s1311(double arad,double aepsge,double amax,int *jstat)
#else
double s1311(arad,aepsge,amax,jstat)
	     double arad;
	     double aepsge;
	     double amax;
	     int    *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To make the step lenth in the iteration procedure based
*              on radius of curvature and an absolute tolerance
*
*
*
* INPUT      : arad   - Radius of curvature
*              aepsge - Absolute tolerance describing the deviation between
*                       the circle of curvature and an Hermite approximation
*                       to the circle
*              amax   - Upper bound of absolute value of coordinates.
*                       If amax = 0.0 is ignored and amax < 0.0 causes error.
*
* OUTPUT     : s1311  - Actual step length to be employed
*              jstat  - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : First a formula giving the angle of the circle segment
*              to be approximated by an Hermit curve within a tolerance
*              is used. Then the arc length of this circular piece is
*              calculated. We make sure that the step length is maximum
*              half the radius of curvature.
*              Two special cases might occure:
*               - Radius of curvature 0
*               - Radius of curvature infinit
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
  int    kpos=1;            /* Position of error                 */
  double tstep;             /* Preliminary value for step length */
  double t1sixth;           /* The value of 1/6                  */
  double talfa;             /* Angle                             */

  if (amax < DZERO) goto err177;

  if (aepsge < DZERO) goto err120;

  if (arad > DZERO)
    {
      t1sixth = (double)1.0/(double)6.0;
      /*  Estimat the opening angle of the segments based on the error
       *   formula. */
      talfa = PI*pow(aepsge/arad,t1sixth)/((double)0.4879);

      /*  Estimate step length equal to curve length of this circular arc,
       *   We limit the step length to half the radius of curvature  */

      tstep = MIN(fabs(talfa*arad),fabs(arad/(double)2.0));
    }
  else if (DEQUAL(arad,DZERO))
    {
      /*  Radius of curvature is zero */
      tstep = (double)100.0*aepsge;
    }

  else
    {
      /*  Infinit radius of curvatur  */
      tstep = amax;
    }

  if ( amax > DZERO && amax < tstep )
    tstep = MAX(amax,aepsge);

  tstep = MAX(tstep,aepsge);

  *jstat = 0;
  goto out;

/* Negative tolerance */

err120: *jstat = -120;
        s6err("s1311",*jstat,kpos);
goto out;

/* Maximal step length zero are less than geometry tolerance */

err177: *jstat = -177;
        s6err("s1311",*jstat,kpos);
goto out;

out:
return(tstep);
}
