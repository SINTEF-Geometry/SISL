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
 * $Id: s2501.c,v 1.1 1995-01-18 09:48:53 pfu Exp $
 *
 */


#define S2501

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2501(SISLSurf *surf, int der, double parvalue[], double derive[],
      double normal[], int *leftknot1, int *leftknot2, double *gaussian,
      int *istat)
#else
 void s2501(surf, der, parvalue, derive, normal, leftknot1, leftknot2, gaussian,
	    istat)
      SISLSurf *surf;
      int der;
      double parvalue[];
      double derive[],
      double normal[],
      int *leftknot1;
      int *leftknot2;
      double *gaussian;
      int *istat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the Gaussian K(u,v) of a Surface for given
*                  values (u,v). This is a lower level routine, used
*                  for evaluation of many K(u,v)'s.
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          der      - Not used.
*      parvalue     - Parameter-value at which to evaluate. Dimension of
*                     parvalue is 2.
*       derive      - Array containing derivatives from routine s1421().
*                     Size = idim*6.
*       normal      - Array containing the normal from routine s1421().
*                     Size = 3.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     first parameter direction where parvalue[0] is found,
*                     that is:
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the routine.
*
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     second parameter direction where parvalue[1] is found,
*                     that is:
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the routine.
*
*  OUTPUT       :
*     gaussian      - Gaussian of the surface in (u,v) =
*                     (parvalue[0],parvalue[1]).
*        istat      - Status messages
*
*                         > 0 : Warning
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The Gaussian is given by
*
*                      K(x,y) = (hxx*hyy-hxy^2)/((1+hx^2+hy^2)^2),
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      K(u,v) = (eg-f*f)/(EG-F*F),
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Gaussian
*                    K(u,v).
*               (ii) If the surface is closed to degenerate, the Gaussian
*                    K(u,v) can be numerical unstable.
*              (iii) The surface should be C2, since the Gaussian is calculated
*                    from the second derivatives. But since the routine is using
*                    right derivatives, the Gaussian will be correct (provided
*                    that the surface is not degenerate).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.  The routine return istat < 0.
*
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.            Date: 1995-1
*****************************************************************************
*/
{
  double a,b;            /* Dummy variables.                                */
  double hx,hy,
    hxx,hyy,hxy;        /* The derivatives of the 1D surface, h(x,y).      */
  double E,F,G;          /* The coefficents of the first fundamental form,
			    that is, E = <Xu,Xu>, F = <Xu,Xv>  and
			    G = <Xv,Xv>.                                    */
  double e,f,g;          /* The coefficents of the second fundamental form,
			    that is, e = <N,Xuu>, f = <N,Xuv> and
			    g = <N,Xvv>.                                   */



  if (surf->idim == 1) /* 1D surface */
  {
    hx  = derive[1];
    hy  = derive[2];
    hxx = derive[3];
    hxy = derive[4];
    hyy = derive[5];

    a = (1+hx*hx+hy*hy);
    *gaussian = (hxx*hyy-hxy*hxy)/(a*a);
  }

  else if (surf->idim == 2) /* 2D surface */
  {
    /* The surface lies in a plane => K(u,v) = 0 */

    *gaussian = 0.0;
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

    /* Compute gaussian = (e*g-f*f)/(E*G-F*F). */
    *gaussian = (e*g-f*f)/(b*b);

  }
  else /* When surf->idim != 1,2 or 3 */
  {
    goto err105;
  }




  /* Successful computations  */

  *istat = 0;
  return;


   /* Error in input, surf->idim != 1,2 or 3 */
err105:
  *istat = -105;
  s6err("s2501",*istat,0);
  return;

}
