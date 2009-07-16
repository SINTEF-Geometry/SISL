/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1995 by                                                     */
/*     SINTEF, Oslo, Norway                                                  */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"

/*
 *
 * $Id: s2542.c,v 1.3 2001-03-19 15:59:00 afr Exp $
 *
 */


#define S2542

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
   s2542(SISLSurf *surf, int ider, int iside1, int iside2, double parvalue[],
      int *leftknot1,int *leftknot2, double *k1, double *k2,  
      double d1[], double d2[],int *jstat)
#else
 void s2542(surf, ider, iside1, iside2,  parvalue, leftknot1, leftknot2, 
	    k1, k2, d1, d2, jstat)
      SISLSurf *surf;
      int    ider;
      int    iside1;
      int    iside2;
      double parvalue[];
      int *leftknot1;
      int *leftknot2;
      double *k1;
      double *k2;
      double d1[];
      double d2[];
      int *jstat;
#endif
/*
***************************************************************************
*
***************************************************************************
*  PURPOSE      :  To compute principal curvature (k1,k2) with corresponding 
*                  principal directions (d1,d2) of a surface for 
*                  given values (u,v) = (parvalue[0],parvalue[1]), where:
*
*                          etl[leftknot1] <= parvalue[0] < etl[leftknot1+1],
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*
*  INPUT        :
*          surf     - Pointer to the surface to evaluate.
*          ider     - Number of derivatives to calculate.
*                     Only implemented for ider=0.
*                       < 0 : No derivative calculated.
*                       = 0 : Principal curvature calculated.
*                       = 1 : Principal curvature and its first derivative 
*                             calculated.
*                       etc.
*          iside1   - Indicator telling if the principal curvature in the first
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate principal curvature from the left hand 
*                             side
*                        >= 0 calculate principal curvature from the right hand 
*                             side.
*          iside2   - Indicator telling if the principal curvature in the second
*                     parameter direction is to be calculated from the
*                     left or from the right:
*                        <  0 calculate principal curvature from the left hand 
*                             side
*                        >= 0 calculate principal curvature from the right hand 
*                             side.
*      parvalue     - Parameter-value at which to evaluate. Dimension of
*                     parvalue is 2.
*
*  INPUT/OUTPUT :
*     leftknot1     - Pointer to the interval in the knot vektor in the 
*                     first parameter direction where parvalue[0] is found,
*                     that is: 
*                          etl[leftknot1] <= parvalue[0] < etl[leftknot1+1].
*                     leftknot1 should be set equal to zero at the first call
*                     to the rutine.
*          
*     leftknot1     - Pointer to the interval in the knot vektor in the 
*                     second parameter direction where parvalue[1] is found,
*                     that is: 
*                          et2[leftknot2] <= parvalue[1] < et2[leftknot2+1].
*                     leftknot2 should be set equal to zero at the first call
*                     to the rutine.
*
*  OUTPUT       :
*         k1        - Max. principal curvature.
*         k2        - Min. principal curvature.
*         d1        - Max. direction of the principal curvature k1, given 
*                     in local coordiantes (with regard to Xu,Xv).
*                     Dimension = 2.
*         d2        - Min. direction of the principal curvature k2, given 
*                     in local coordiantes (with regard to Xu,Xv).
*                     Dimension = 2.
*        jstat      - Status messages
*                         = 2 : Surface is degenerate at the point, that is,
*                               the surface is not regular at this point.
*                         = 1 : Surface is close to degenerate at the point.
*                               Angle between tangents is less than the angular
*                               tolerance.
*                         = 0 : Ok.
*                         < 0 : Error.
*
*
*  METHOD        :  The princpal curvatures -k1 and -k2 are eigenvalues of
*                   dN, thus it turn out that we have to solve a 
*                   eigenvalue/eigenvector problem, see references.
*
*  REFERENCES   :  Differential Geometry of Curves and Surfaces,
*                    (Manfredo P. Do Carmo, Prentice Hall,
*                      ISBN: 0-13-212589-7).
*
*                  Elementary Linear Algebra 5e
*                    (Howard Anton, Wiley, ISBN:0-471-84819-0)
*-
*  CALLS        :  s1422() and s2543().
*
*  LIMITATIONS  :
*                (i) If the surface is degenrated (not regular) at the point
*                    (u,v), it makes now sence to speak about curvature.
*                    The routine return istat == 2.
*               (ii) If the surface is close to degenrate, the resuts
*                    can be numerical unstable. The routine return 
*                    istat == 1.
*              (iii) The dimension of the space in which the surface lies must 
*                    be 1,2 or 3.
*  
*
* WRITTEN BY    :  Geir Westgaard,  SINTEF, Oslo, Norway.         Date: 1995-1
******************************************************************************
*/
{
  int kwarn = 0;         /* Local staus variable(warning).                  */
  int der = 0;           /*  (dummy) */
  double derive[18];     /* Array containing the computed derivatives.      */
  double normal[3];      /* Array containing the computed normalvektor.     */


  if (ider != 0) goto err178;


  if (surf == SISL_NULL)  goto err150;
  else
  {
    /* Compute derivates and normal. */

    s1422(surf,2,iside1,iside2,parvalue,leftknot1,leftknot2,derive,normal,
	  jstat);
    if (*jstat > 0) kwarn = *jstat;

    if (*jstat < 0) /* Error in lower level routine. */
    {
      goto error;
    }
    else if (*jstat != 2) /* The surface is not degenerate */
    {
      s2543(surf, der, derive, normal, k1, k2, d1, d2, jstat);

      if (*jstat < 0)
	goto error;
    }
    else if (*jstat == 2) /* The surface is degenerated. */
    {
      *k1 = 0.0;
      *k2 = 0.0;
      d1[0] = 1.0;
      d1[1] = 0.0;
      d2[0] = 0.0;
      d2[1] = 1.0;
      
      goto war002;
    }

  }
  
  /* Successful computations  */
  
  *jstat = kwarn;
  goto out;
    
  /* The surface is degenerated at (u,v) */
  war002:
     *jstat = 2;
  goto out;
  
  /* Error. Input (surface) pointer is SISL_NULL. */
  err150:
     *jstat = -150;
  s6err("s2542", *jstat, 0);
  goto out;
  
  /* Illegal derivative requested. */
  err178:
     *jstat = -178;
  s6err("s2542",*jstat,0);
  goto out;
  /* Error in lower level routine.  */
  
  error:
     s6err("s2542",*jstat,0);
  goto out;
  
  
  out:
     
     return;

}


