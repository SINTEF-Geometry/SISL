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
 * $Id: s1332.c,v 1.3 2001-03-19 15:58:45 afr Exp $
 *
 */


#define S1332

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
s1332(SISLCurve *pc1,SISLCurve *pc2,double aepsge,double ep[],SISLSurf **rs,int *jstat)
#else
void s1332(pc1,pc2,aepsge,ep,rs,jstat)
     SISLCurve  *pc1;
     SISLCurve  *pc2;
     double aepsge;
     double ep[];
     SISLSurf   **rs;
     int *jstat;
#endif
/*
*********************************************************************
*
* PURPOSE    : To create a swept B-spline surface by making
*              the tensor-product of two B-spline curves.
*
* INPUT      : pc1    - Pointer curve 1.
*              pc2    - Pointer curve 2.
*              aepsge - Maximal deviation allowed between true swept
*                       surface and generated surface.
*              ep     - SISLPoint near the curve to be swept. If the point
*                       lies on curve 2, then curve 2 is swept along curve 1
*                       with the point as contact point. If the point lies
*                       on curve 1, then curve 1 is swept along curve 2 with
*                       the point as contact point. If the point is not lying
*                       on any of the curve, then the surface will not
*                       interpolate any of the curves.
*
*
* OUTPUT     : jstat  - status messages
*                        > 0      : warning
*                        = 0      : ok
*                        < 0      : error
*              rs     - Pointer to the surface produced.
*
* METHOD     : Mathematically, the surface is expressed as:
*
*                       in1 in2
*              P(u,v) = SUM SUM ((p1(i) + p2(j)-ep)*w1(i)*w2(j)*B(i,ik1)*B(j,ik2)
*                       i=1 j=1
*-
* CALLS      : s1707,s6err.
*
* WRITTEN BY : A. M. Ytrehus  SI,  Oslo, Norway.  Sep. 1988
* REVISED BY : Johannes Kaasa SI,  Sep. 1991 (Introduced NURBS)
* Revised by : Paal Fugelli, SINTEF, Oslo, Norway, Oct. 1994. Removed memory leaks.
*********************************************************************
*/
{
  double *sknot1 = SISL_NULL;  /* Pointer to knot-vector of curve 1.         */
  double *scoef1 = SISL_NULL;  /* Pointer to vertices of curve 1.            */
  double *rcoef1 = SISL_NULL;  /* Pointer to rational vertices of curve 1.   */
  int kord1;              /* Order of curve 1.                          */
  int kn1;                /* Number of vertices of curve 1.             */
  double *sknot2 = SISL_NULL;  /* Pointer to knot-vector of curve 2.         */
  double *scoef2 = SISL_NULL;  /* Pointer to vertices of curve 2.            */
  double *rcoef2 = SISL_NULL;  /* Pointer to rational vertices of curve 2.   */
  int kord2;              /* Order of curve 2.                          */
  int kn2;                /* Number of vertices of curve 2.             */
  int kdim;               /* Dimension of the space in which the
			     curves lies.                               */
  int rdim;               /* Dimension of rational space.               */
  double *scoef = SISL_NULL;   /* Pointer to vertex array for surface        */
  double *sp,*spc1,*spc2; /* Pointers to scoef, scoef1 snd scoef2.      */
  double *spnt;           /* Pointer to epoint.                         */
  int ki,kj,kp;           /* Loop controllers.                          */
  int kcopy,kind;         /* Variables for use in newSurf.              */
  int kstat = 0;          /* Status variable                            */
  int kpos = 0;           /* Position of error                          */

  double *weight1 = SISL_NULL; /* Rational weights in first direction        */
  double *weight2 = SISL_NULL; /* Rational weights in second direction       */
  double weight;          /* Tensor product weights                     */


  *rs = SISL_NULL;

  /* The curves must have the same dimension      */
  if (pc1 -> idim != pc2 -> idim) goto err106;

  /* Make local pointers to description of curves */

  if (!pc1) goto err150;
  if (!pc2) goto err150;

  s1707(pc1,&kstat);
  if (kstat<0) goto error;

  s1707(pc2,&kstat);
  if (kstat<0) goto error;

  sknot1 = pc1 -> et;
  scoef1 = pc1 -> ecoef;
  rcoef1 = pc1 -> rcoef;
  kn1 = pc1 -> in;
  kord1 = pc1 -> ik;
  kdim = pc1 -> idim;

  sknot2 = pc2 -> et;
  scoef2 = pc2 -> ecoef;
  rcoef2 = pc2 -> rcoef;
  kn2 = pc2 -> in;
  kord2 = pc2 -> ik;

  /* Allocate vertex-array for the surface. */

  rdim = kdim + 1;
  if (pc1->ikind == 2 || pc1->ikind == 4
                      || pc2->ikind == 2 || pc2->ikind == 4)
    scoef = newarray(kn1*kn2*rdim,DOUBLE);
  else
    scoef = newarray(kn1*kn2*kdim,DOUBLE);
  if (!scoef) goto err101;

  /* Allocate and initiate rational weights. */

  weight1 = newarray(kn1,DOUBLE);
  if (!weight1) goto err101;
  if (pc1->ikind == 2 || pc1->ikind == 4)
    for (ki=0; ki<kn1; ki++)
      weight1[ki] = rcoef1[(ki + 1)*rdim - 1];
  else
    for (ki=0; ki<kn1; ki++)
      weight1[ki] = 1.;
  weight2 = newarray(kn2,DOUBLE);
  if (!weight2) goto err101;
  if (pc2->ikind == 2 || pc2->ikind == 4)
    for (ki=0; ki<kn2; ki++)
      weight2[ki] = rcoef2[(ki + 1)*rdim - 1];
  else
    for (ki=0; ki<kn2; ki++)
      weight2[ki] = 1.;

  /* Compute the vertices of the surface. */

  sp = scoef;
  for (kj=0; kj<kn2; kj++)
    {
      for (ki=0; ki<kn1; ki++)
        {
	  spc1 = scoef1 + ki*kdim;
	  spc2 = scoef2 + kj*kdim;
	  spnt = ep;
          weight = weight1[ki]*weight2[kj];

	  for (kp=0; kp<kdim; kp++)
            {
	      *sp = (*spc1 + *spc2 - *spnt)*weight;
	      sp++;
	      spc1++;
	      spc2++;
	      spnt++;
            }
          if (pc1->ikind == 2 || pc1->ikind == 4
                      || pc2->ikind == 2 || pc2->ikind == 4)
            {
              *sp = weight;
              sp++;
            }
        }
    }

  /* Create the surface */

  kcopy = 1;
  if (pc1->ikind == 2 || pc1->ikind == 4
                      || pc2->ikind == 2 || pc2->ikind == 4)
    kind = 2;
  else
    kind = 1;
  *rs =  newSurf(kn1,kn2,kord1,kord2,sknot1,sknot2,scoef,kind,kdim,kcopy);

  *jstat = 0;
  goto out;

  /* Empty curve. */
 err150: *jstat = -150;
  s6err("s1332",*jstat,kpos);
  goto out;

  /* Different dimension for curves   */
 err106: *jstat = -106;
  s6err("s1332",*jstat,kpos);
  goto out;

  /* Error in space allocation. */
 err101: *jstat = -101;
  s6err("s1332",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */
 error:  *jstat = kstat;
  s6err("s1332",*jstat,kpos);
  goto out;
 out:

  /* Free local memory */
  if (scoef) freearray(scoef);
  if (weight1) freearray(weight1);
  if (weight2) freearray(weight2);

  return;
}
