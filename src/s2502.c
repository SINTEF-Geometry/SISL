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
 * $Id: s2502.c,v 1.2 1995-01-18 13:14:49 pfu Exp $
 *
 */


#define S2502

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
   s2502(SISLSurf *surf, double parvalue[], int *leftknot1,
      int *leftknot2, double *meancurvature, int *istat)
#else
 void s2502(surf, parvalue, leftknot1, leftknot2, meancurvature, istat)
      SISLSurf *surf;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double *meancurvature;
      int *istat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the mean curvature H(u,v) of a surface for given
*                  values (u,v) = (parvalue[0],parvalue[1]), where:
*
*                          etl[leftknot1] <= parvalue[0] < etl[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*      parvalue     - Parameter-value at which to evaluate. Dimension of
*                     parvalue is 2.
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
*    meancurvature  - Mean curvature of the surface in (u,v) =
*                     (parvalue[0],parvalue[1]).
*        istat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD        :  The mean curvature is given by
*
*                      H(x,y) = 0.5((1+hx^2)hyy-2hxhyhxy+(1+hy^2)hxx)/
*                                      ((1+hx^2+hy^2)^3/2),
*
*                  if the surface (h(x,y)) is 1D, and
*
*                      H(u,v) = 0.5(eG-2fF+gE)/(EG-F*F),
*
*                  if the surface (X(u,v)) is 3D. The variables E,F,G,e,f and g
*                  are the coefficients of the first and second fundamental form.
*                  They are given by: e = <N,Xuu>, f = <N,Xuv>, g = <N,Xvv>,
*                  E = <Xu,Xu>, F = <Xu,Xv> and G = <Xv,Xv>. The routine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1421() and s2503().
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the mean curvature
*                    H(u,v). The routine returns istat = 2.
*               (ii) If the surface is closed to degenerate, the mean curvature
*                    H(u,v) can be numerical unstable. The routine returns
*                    istat = 1.
*              (iii) The surface should be C2, since the mean c. is calculated
*                    from the second derivatives. But since the routine is using
*                    right derivatives, the mean c. will be correct (provided
*                    that the surface is not degenerate).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.  The routine returns istat < 0.
*
*
* WRITTEN BY :  Geir Westgaard, SINTEF, Oslo, Norway.             Date: 1995-1
******************************************************************************
*/
{
  int kistat = 0;     /* Local staus variable.                           */
  int der = 0;        /*  (dummy) */
  double derive[18];  /* Array containing the computed derivatives.      */
  double normal[3];   /* Array containing the computed normalvektor.     */




  if (surf == NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

    s1421(surf,2,parvalue,leftknot1,leftknot2,derive,normal,&kistat);

    if (kistat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (kistat != 2) /* The surface is not degenerate */
    {
      s2503(surf, der, parvalue, derive, normal, leftknot1, leftknot2,
	    meancurvature, &kistat);

      if (kistat < 0)
	goto error;
    }
    else if (kistat == 2) /* The surface is degenerated. */
    {
      *meancurvature = 0.0;
      goto war002;
    }

  }


  /* Successful computations  */

  *istat = kistat;
  goto out;



   /* The surface is degenerated at (u,v) */
war002:
  *istat = 2;
  goto out;

  /* Error. Input (surface) pointer is NULL. */
err150:
  *istat = -150;
  s6err("s2502", *istat, 0);
  goto out;

   /* Error in lower level routine.  */
error:
  *istat = kistat;
  s6err("s2502",*istat,0);
  goto out;

out:

  return;
}
