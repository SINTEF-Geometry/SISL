/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1995 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: s2505.c,v 1.4 1995-06-28 11:17:48 jka Exp $
 *
 */


#define S2505

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2505(SISLSurf *surf, int der, double derive[], double normal[],
      double *absCurvature, int *jstat)
#else
 void s2505(surf, der, derive, normal, absCurvature, jstat)
      SISLSurf *surf;
      int der;
      double derive[],
      double normal[],
      double *absCurvature;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the absolute curvature K(u,v) of a Surface
*                  for given values (u,v). This is a lower level routine,
*                  used for evaluation of many K(u,v)'s.
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          der      - Not used.
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  INPUT/OUTPUT :
*
*  OUTPUT       :
*    absCurvature   - Absolute curvature of the surface in (u,v) =
*        jstat      - Staus messages
*                         = 2 : Surface is degenrate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is closed to degenrate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : error.
*
*  METHOD        :  The absolute curvature is given by
*
*                      A(x,y) = |k1| + |k2|,
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      A(u,v) = |k1| + |k2|,
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficents of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The rutine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1421() and s1424().
*
*  LIMITATIONS  :
*                (i) If the surface is degenrated (not regular) at the point
*                    (u,v), it makes now sence to speak about the absolute c.
*                    A(u,v). The routine return jstat == 2.
*               (ii) If the surface is closed to degenrate, the absolute c.
*                    A(u,v) can be numerical unstable. The routine return
*                    jstat == 1.
*              (iii) The surface should be C2, since the absolute c. is calculated
*                    from the second derivativs. But since the rutine is using
*                    right derivatives, the absolute c. will be correct (provided
*                    that the surface is not degenerate).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.  The routine return jstat < 0.
*
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.            Date: 1995-1
*****************************************************************************
*/
{
  double a,b;             /* Temporary variables.                                */
  double hx,hy,
    hxx,hyy,hxy;         /* The derivatives of the 1D surface, h(x,y).      */
  double E,F,G;           /* The coefficents of the first fundamental form,
			     that is, E = <Xu,Xu>, F = <Xu,Xv>  and
			     G = <Xv,Xv>.                                    */
  double e,f,g;          /* The coefficents of the second fundamental form,
			    that is, e = <N,Xuu>, f = <N,Xuv> and
			    g = <N,Xvv>.                                     */
  double gc;             /* Gaussian curvature.                              */
  double mc;             /* Mean curvature.                                  */



  if (surf->idim == 1) /* 1D surface */
  {
    hx  = derive[1];
    hy  = derive[2];
    hxx = derive[3];
    hxy = derive[4];
    hyy = derive[5];


    a = (1+hx*hx+hy*hy);
    gc = (hxx*hyy-hxy*hxy)/(a*a);

    b = sqrt(a*a*a);

    a = (1.0 + hx*hx)*hyy - 2.0*hx*hy*hxy + (1.0 + hy*hy)*hxx;

    mc = 0.5*a/b;

    *absCurvature = fabs(mc + sqrt(mc*mc - gc)) +
       fabs(mc - sqrt(mc*mc - gc));
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => A(u,v) = 0 */

    *absCurvature = 0.0;
  }
  else if (surf->idim == 3) /* 3D surface */
  {
    /* E = <Xu,Xu> */
    E = derive[3]*derive[3]+derive[4]*derive[4]+derive[5]*derive[5];

    /* F = <Xu,Xv> */
    F = derive[3]*derive[6]+derive[4]*derive[7]+derive[5]*derive[8];

    /* G = <Xv,Xv> */
    G = derive[6]*derive[6]+derive[7]*derive[7]+derive[8]*derive[8];

    /* b = EG + F^2. */
    b = E*G-F*F;

    /* e = <N,Xuu> (/ sqrt(E*G-F*F)) */
    e = normal[0]*derive[9]+normal[1]*derive[10]+normal[2]*derive[11];

    /* f = <N,Xuv> (/ sqrt(E*G-F*F)) */
    f = normal[0]*derive[12]+normal[1]*derive[13]+normal[2]*derive[14];

    /* g = <N,Xvv> (/ sqrt(E*G-F*F)) */
    g = normal[0]*derive[15]+normal[1]*derive[16]+normal[2]*derive[17];

    /* Compute absolute curvature = |k1| + |k2|. */

    gc = (e*g-f*f)/(b*b);

    a = 0.5*(e*G - 2.0*f*F + g*E);
    b = b*sqrt(b);

    mc = a/b;

    *absCurvature = fabs(mc + sqrt(mc*mc - gc)) +
       fabs(mc - sqrt(mc*mc - gc));
  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }




  /* Successful computations  */

  *jstat = 0;
  goto out;


  /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *jstat = -105;
  s6err("s2505",*jstat,0);
  goto out;

out:

  return;

}
