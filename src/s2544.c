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
 * $Id: s2544.c,v 1.2 1995-09-22 13:23:26 jka Exp $
 *
 */


#define S2544

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s2544(SISLSurf *surf, int ider, int iside1, int iside2, double parvalue[],
      int *leftknot1,int *leftknot2, double *norcurv, int *jstat)
#else
 void s2544(surf, ider, iside1, iside2, parvalue, leftknot1, leftknot2, norcurv,
	    jstat)
      SISLSurf *surf;
      int    ider;
      int    iside1;
      int    iside2;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double *norcurv;
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute the Normal curvature of a Surface for given
*                  values (u,v) = (parvalue[0],parvalue[1]), in the
*                  direction (parvalue[2],parvalue[3])
*                  where:
*
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                  
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0.
*                       < 0 : No derivative calculated.
*                       = 0 : Position calculated.
*                       = 1 : Position and first derivative calculated.
*                       etc.
*          iside1   - Indicator telling if the derivatives in the first
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*          iside2   - Indicator telling if the derivatives in the second
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate derivative from the left hand side
*                        >= 0 calculate derivative from the right hand side.
*      parvalue     - Parameter-value at which to evaluate pluss the direction
*                     Dimension of parvalue is 4.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vector in the
*                     first parameter direction where parvalue[0] is found,
*                     that is:
*                          et1[leftknot1] <= parvalue[0] < et1[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the routine.
*
*     leftknot2     - Pointer to the interval in the knot vector in the
*                     second parameter direction where parvalue[1] is found,
*                     that is:
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the routine.
*
*  OUTPUT       :
*     norcurv      - Normal curvature of the surface in (u,v) =(parvalue[0],parvalue[1])
*                    in the direction (parvalue[2],parvalue[3]).
*        jstat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*  METHOD       :  The Normal curvature is given by
*
*                      <Xuu,N>d1*d1 + 2<Xuv,N>d1*d2 + <Xvv,N>d2*d2
*                      --------------------------------------------
*                      <Xu,Xu>d1*d1 + 2<Xu,Xv>d1*d2 + <Xv,Xv>d2*d2
*                  
*                  
*		   Where (d1,d2) is the normalised parameter direction
*                  (parvalue[2],parvalue[3]) and the other numbers are
*                  the factors for the first and second fundamental forms.
*
*                  The routine will
*                  test if the surface is degenerate (not regular) or close to
*                  degenerate. 
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*-
*  CALLS        :  s1422() and s2513().
*
*  LIMITATIONS  :
*                (i) If the surface is degenerated (not regular) at the point
*                    (u,v), it makes no sense to speak about the Normal Curvature.
*                    The routine returns jstat = 2.
*               (ii) If the surface is close to degenerate, the calculations
*                    can be numerical instable. The routine returns
*                    jstat = 1.
*              (iii) If the surface is Cr the curvature calculated is C(r-2).
*               (iv) The dimension of the space in which the surface lies must
*                    be 1,2 or 3.
*
*
* WRITTEN BY :    Ulf J Krystad, SINTEF, Oslo, Norway.            Date: 1995-1
* REVISED BY :    Johannes Kaasa, SINTEF, Oslo, Norway.           Date: 1995-9
******************************************************************************
*/
{
  int kwarn = 0;      	 /* Local staus variable(warning).                  */
  double derive[18];     /* Array containing the computed derivatives.      */
  double normal[3];      /* Array containing the computed normalvektor.     */
  double val[6];         /* Array containing the computed values in s2511   */
  double d1, d2, dl;     /* The normalised parameter direction
			    (parvalue[2],parvalue[3]) and their length      */
  double temp1, temp2;   /* Temporary values				    */
  /* ______________________________________________________________________ */
  

  if (ider != 0) goto err178;
  
  dl = sqrt(parvalue[2]*parvalue[2] + parvalue[3]*parvalue[3]);
  if (dl < REL_PAR_RES) goto err174;
  d1 = parvalue[2]/dl;
  d2 = parvalue[3]/dl;


  if (surf == NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

    s1422(surf, 2, iside1, iside2, parvalue, leftknot1, leftknot2, derive,
	  normal, jstat);
    if (*jstat > 0) kwarn = *jstat;

    if (*jstat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (*jstat != 2) /* The surface is not degenerate */
    {
       /* Find factors in fundamental form */
       s2513(surf, 0, 2, 1, derive, normal, val, jstat);
       
       if (*jstat < 0)
	  goto error;
       
       temp1    = val[0]*d1*d1 + 2*val[1]*d1*d2 + val[2]*d2*d2;
       if (temp1 < REL_PAR_RES) goto err174;
       temp2    = val[3]*d1*d1 + 2*val[4]*d1*d2 + val[5]*d2*d2;
       
       *norcurv = temp2/temp1;

       
    }
    else if (*jstat == 2) /* The surface is degenerated. */
    {
      *norcurv = 0.0;
      goto war002;
    }

  }


  /* Successful computations  */

  *jstat = kwarn;
  goto out;


  /* ____________________________________________________________________ */
  /*                           ERROR EXIT				  */
  /* ____________________________________________________________________ */
  
   /* The surface is degenerated at (u,v) */
war002:
  *jstat = 2;
  goto out;

  /* Error. Input (surface) pointer is NULL. */
err150:
  *jstat = -150;
  s6err("s2544", *jstat, 0);
  goto out;

  /* Degenerate condition. */
err174:
  *jstat = -174;
  s6err("s2544",*jstat,0);
  goto out;
  
  /* Illegal derivative requested. */
err178:
  *jstat = -178;
  s6err("s2544",*jstat,0);
  goto out;
  /* Error in lower level routine.  */

error:
  s6err("s2544",*jstat,0);
  goto out;


  /* ____________________________________________________________________ */
  /*                        THE ONE AND ONLY EXIT			  */
  /* ____________________________________________________________________ */
out:

  return;

}
