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
 * $Id: s1307.c,v 1.2 2001-03-19 15:58:43 afr Exp $
 *
 */
#define S1307

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1307(double ep[],int idim,double egeo[],int *jstat)
#else
void s1307(ep,idim,egeo,jstat)
     double ep[];
     int    idim;
     double egeo[];
     int    *jstat;
#endif
/*
*********************************************************************
*                                                                   
* PURPOSE    : To calculate the unit tangent,curvature and radius of
*              curvature of a curve at a point.
*
* INPUT      : ep     - Position, first and second derivative of the
*                       curve with respect to some parametrization 
*                       at the point. (3*idim doubles)
*              idim   - Dimension of the space the curve lies in
*
*
* OUTPUT     : 
*              jstat  - status messages  
*                         = 1      : Curvature radius infinit
*                         = 0      : ok, curvature radius
*                         < 0      : error
*              egeo   - 3-D geometry description of the intersection. The
*                       array contains: position, unit tangent, curvature
*                       and radius of curvature. (A total of 3*idim + 1
*                       doubles). A radius of curvature =-1, indicates
*                       that the radius of curvature is infinit.
*
* METHOD     : We convert the description of the derivatives to an
*              arc length parametrization. In this parametrization
*              the second derivative is the same as the curvature vector.
*              The radius of curvature is the invers of the length of this
*              curvature vector.
*
* REFERENCES : 
*-
* CALLS      : s6scp,s6norm,s6length
*              
*
* WRITTEN BY : Tor Dokken, SI, Oslo , Norway, 3 July 1988
* Revised by : Tor Dokken, SI, Oslo , Norway, March 1989
*              Corrected use of maximal radius of curvature.
*
*********************************************************************
*/
{
  int k2dim=2*idim;   /* The dimension *2, Start of double derivative*/
  int kstat;          /* Local status variable                       */
  int ki,kj;          /* Variables in loop                           */
  double tlength;     /* Length of first derivative vector           */
  double tdum;        /* Dummy variable                              */
  
  /* Let c = c(w) be a parameterized curve.
   *  The curvature vector is defined as the derivative of the unit tangent
   *  vector with respect to the arc length a. If we don't have an arclength
   *  parametrization then this parametrization can be written as a function
   *  of the arc length w = w(a). By using the kernel rule for differentiation
   *  we get:
   *
   *         d            d       dw   d    c'(w)    dw   d    c'(w)      da
   *  k(a) = -- T(w(a)) = -- T(w) -- = -- ---------- -- = -- ---------- / --
   *         da           dw      da   dw sqrt(c'c') da   dw sqrt(c'c')   dw
   *
   *
   *         d       c'(w)                c"        c' (c'c'')
   *         -- ----------------- =   ---------- - ------------- 
   *         dw sqrt(c'(w) c'(w))     sqrt(c'c')   sqrt(c'c')**3
   *
   *
   *
   *         da
   *         -- = sqrt(c'c')
   *         dw 
   */
  
  /* Copy position */
  
  memcopy(egeo,ep,idim,DOUBLE);
  
  /* First we normalize the tangent vector */
  
  tlength = s6norm(ep+idim,idim,egeo+idim,&kstat);
  
  if (DEQUAL(tlength,(double)0.0)) goto war101;
  
  /* Make curvature vector */
  
  tdum = s6scpr(ep+k2dim,egeo+idim,idim)/tlength;
  
  for (ki=idim,kj=k2dim;ki<k2dim;ki++,kj++)
    {
      egeo[kj] = (ep[kj]/tlength - egeo[ki]*tdum)/tlength;
    }
  
  /* Make radius of curvature */
  
  tdum = s6length(egeo+k2dim,idim,&kstat);
  
  if (tdum!=DZERO && ((double)1.0/tdum) > MAXIMAL_RADIUS_OF_CURVATURE) 
    goto war101;
  
  if (DNEQUAL(tdum,(double)0.0))
    {
      egeo[3*idim] = (double)1.0/tdum;
    }
  else
    {
      goto war101;
    }
  
  /* Everyting is ok */
  
  *jstat = 0;
  goto out;
  
  /* Infinit radius of curvature */
  
 war101: *jstat=1;
  egeo[3*idim] = (double)-1.0;
  goto out;
  
 out:
  return;
}
