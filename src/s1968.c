/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the sisl-copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "sisl-copyright.h"


#define S1968

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1968(SISLSurf *ps,double eeps[],int nend[],int iopen1,int iopen2,
	   double edgeps[],int iopt,int itmax,SISLSurf **rs,
	   int *jstat)
#else
void s1968(ps,eeps,nend,iopen1,iopen2,edgeps,iopt,itmax,rs,jstat)
     SISLSurf   *ps;
     double eeps[];
     int    nend[];
     int    iopen1;
     int    iopen2;
     double edgeps[];
     int    iopt;
     int    itmax;
     SISLSurf   **rs;
     int    *jstat;
#endif
/*
********************************************************************
*
* Purpose: To compute a cubic tensor-product spline approximation
*          to a given tensor product spline surface of arbitrary order,
*          with error less than eeps in each of the idim components.
*          Mark that the error in continuity over the start and end of
*          a closed or periodic surface is only guaranteed to be within
*          edgeps.
*
* Input : Ps     - Pointer to surface.
*         Eeps   - Array (length kdim) containing the maximum deviation
*                  which is acceptable in each of the kdim components of
*                  of the surface (except possibly along the edges).
*         Nend   - Array (length 4) containing the no. of derivatives
*                  to be kept fixed along each edge of the surface.
*                  The numbering of the edges is the same as for edeps below.
*                  All the derivatives of order < nend(i)-1 will be kept fixed
*                  along edge no. i. Hence nend(i)=0 indicates that nothing
*                  is to be kpet fixed along edge no. i.
*                  To be kept fixed here means to have error less than edgeps.
*                  In general, it is impossible to remove any knots and
*                  keep an edge completely fixed.
*         iopen1 - Open/closed parameter in first parameter direction.
*                      =  1 : Produce open surface.
*                      =  0 : Produce closed, non-periodic surface if possible.
*                      = -1 : Produce closed, periodic surface if possible.
*         iopen2 - Open/closed parameter in second parameter direction.
*                      =  1 : Produce open surface.
*                      =  0 : Produce closed, non-periodic surface if possible.
*                      = -1 : Produce closed, periodic surface if possible.
*         Edgeps - Array (length kdim*4) containing the max. deviation
*                  which is acceptable along the edges of the surfaces.
*                  Edgeps(1,i):edgeps(kdim,i) gives the tolerance along
*                  the edge corresponding to the i-th parameter having
*                  its min. or max value.
*                  i=1 : min value of first parameter.
*                  i=2 : max value of first parameter.
*                  i=3 : min value of second parameter.
*                  i=4 : max value of second parameter.
*                  Edgeps(kp,i) will only have any significance if nend(i)>0.
*         Iopt   - Flag indicationg the order in which tha data-reduction
*                  is to be performed.
*                   = 1 : Remove knots in parameter-direction 1 only.
*                   = 2 : Remove knots in parameter-direction 2 only.
*                   = 3 : Remove knots in parameter-direction 1 and
*                         and then in parameter-direction 2.
*                   = 4 : Remove knots in parameter-direction 2 and
*                         and then in parameter-direction 1.
*         Itmax  - Max. no. of iterations.
*
* Output:
*         Jstat  - Output status:
*                   < 0 : Error.
*                   = 0 : Ok.
*                   > 0 : Warning:
*         Rs     - Pointer to surface.
*
* Method:
*        First a tensor-product bicubic spline-approximation with error
*        less than half the tolerance is computed. Then the no. of
*        parameters is reduced by the use of data-reduction, also with
*        a tolerance of half the input-tolerance.
*        The initial approximation is deteremined by first sampling the
*        surface and its 10, 01 and 11 derivatives sufficiently often to
*        ensure that a piecewise bicubic hermite interpolation will give
*        an error less than half the tolerance, and then computing this
*        approximation.
*
* Calls: s1387, s1530, s1965, s6err.
*
* Written by: C.R.Birkeland, Si, Oslo, Norway, April 1993.
* The datareduction routine, s1345, is written by: Knut M|rken,  SI.
* Changed by: Per OEyvind, SINTEF, 1994-11.
*             Added input check for invalid Eeps values.
*             Removed 2 array "overlooks", when setting side1 and side2.
**********************************************************************
*/
{
  int i,j,k;                    /* Loop control                        */
  int index = 0;
  int idim=ps->idim;
  int stat = 0;                 /* Error message parameters            */
  int kpos = 0;
  int im1, im2;                 /* Number of vertices                  */
  int leftknot1 = 0;            /* Used in call to s1424: evaluator    */
  int leftknot2 = 0;
  int side1 = 0;                /* Evaluation from what side           */
  int side2 = 0;
  double *error1 = SISL_NULL;        /* For error bounds                    */
  double *error2 = SISL_NULL;
  double *error1472 = SISL_NULL;
  double *newcoeff = SISL_NULL;
  double *ep = SISL_NULL;            /* Used to store points                */
  double *eder10 = SISL_NULL;        /* Used to store corresponding derivatives */
  double *eder01 = SISL_NULL;
  double *eder11 = SISL_NULL;
  double *derive = SISL_NULL;
  double *par1 = SISL_NULL;          /*  Used to store parametrizations      */
  double *par2 = SISL_NULL;
  SISLCurve *ocurve = SISL_NULL;
  SISLSurf *osurf1 = SISL_NULL;
  SISLSurf *osurf2 = SISL_NULL;
  double parvalue[2];


  /* Check input and description of surface.
   * ---------------------------------------
   */

  if (!ps) goto err150;
  if (ps->ik1 < 1 || ps->ik2 < 1 || idim < 1) goto err103;

  /* Get the eeps value with largest absolute value */
  parvalue[0] = fabs(eeps[0]);
  for (i = 1; i < 4; i++)
  {
    parvalue[1] = fabs(eeps[i]);
    if (parvalue[1] > parvalue[0])
      parvalue[0] = parvalue[1];
  }

  /* If no max deviation are large than zero ... */
  if (DEQUAL(parvalue[0], 0.0))
    goto err103;

  /* Allocate error matrices */

  error1 = newarray( idim, DOUBLE);
  error2 = newarray( idim, DOUBLE);
  error1472 = newarray( idim, DOUBLE);
  if( error1 == SISL_NULL || error2 == SISL_NULL || error1472 == SISL_NULL )
    goto err101;

  /* First,  make sure both orders are higher or equal to 4 */

  if (ps->ik1 < 4 || ps->ik2 < 4)
    {
      /* Increase order in at least one direction */

      s1387(ps, MAX( 4, ps->ik1), MAX( 4, ps->ik2),
	    &osurf1, &stat);
      if (stat<0) goto error;
    }
  else
    osurf1 = newSurf(ps->in1, ps->in2, ps->ik1, ps->ik2,
		     ps->et1, ps->et2, ps->ecoef, 1, idim, 1);

  /* Make sure orders in both directions are 4. */

  if (osurf1->ik1 == 4 && osurf1->ik2 == 4)
    {
      /* Orders in both parameterdirections are now 4.
       * No degree reduction necessary.
       * Specify directly error bounds for use
       * in the datareduction routine. */

      for(i=0; i<idim; i++)
	error2[i] = eeps[i];
      osurf2 = newSurf(osurf1->in1, osurf1->in2, 4, 4, osurf1->et1,
		       osurf1->et2, osurf1->ecoef, 1, idim, 1);
    }
  else
    {
      /* The order must be decreased for at least one direction.
       * Create error bounds. */

      for(i=0; i<idim; i++)
	{
	  error1[i] = 0.25 * eeps[i];
	  error2[i] = eeps[i]-error1[i];
	}

      /* Check second direction */

      /* Create curve equivalent and determine sampling points
       * for 2. direction */

      ocurve = newCurve(osurf1->in2, osurf1->ik2, osurf1->et2,
			osurf1->ecoef, 1, idim*osurf1->in1, 1);
      if(ocurve == SISL_NULL) goto err101;

      error1472 = increasearray(error1472, idim*osurf1->in1, DOUBLE);
      if(error1472 == SISL_NULL) goto err101;
      for(i=0; i<idim*osurf1->in1; i++)
	error1472[i] = error1[i % idim];
      s1355(ocurve, error1472, &par2, &im2, &stat);
      if (stat<0) goto error;
      freeCurve(ocurve);
      ocurve = SISL_NULL;

      /* Create curve equivalent and determine sample points
       * for 1. direction after transposing surface coefficients */

      newcoeff = newarray(idim * osurf1->in1 * osurf1->in2, DOUBLE);
      if (newcoeff == SISL_NULL) goto err101;
      s6chpar(osurf1->ecoef, osurf1->in1, osurf1->in2, idim, newcoeff);

      ocurve = newCurve(osurf1->in1, osurf1->ik1, osurf1->et1,
			newcoeff, 1, idim*osurf1->in2, 1);
      if(ocurve == SISL_NULL) goto err101;

      error1472 = increasearray(error1472, idim*osurf1->in2, DOUBLE);
      if(error1472 == SISL_NULL) goto err101;
      for(i=0; i<idim*osurf1->in2; i++)
	error1472[i] = error1[i % idim];
      s1355(ocurve, error1472, &par1, &im1, &stat);
      if (stat<0) goto error;
      freeCurve(ocurve);
      ocurve = SISL_NULL;

      /* Compute points and derivatives on the surface
       * at the points given by calculated par.values */

      derive = newarray( idim * 4, DOUBLE );
      ep     = newarray( idim * im1 * im2, DOUBLE );
      eder10 = newarray( idim * im1 * im2, DOUBLE );
      eder01 = newarray( idim * im1 * im2, DOUBLE );
      eder11 = newarray( idim * im1 * im2, DOUBLE );
      if (ep == SISL_NULL || eder10 == SISL_NULL || eder01 == SISL_NULL ||
	  eder11 == SISL_NULL || derive == SISL_NULL) goto err101;
      index = 0;
      for(j=0; j<im2; j++)
	{
	  parvalue[1] = par2[j];
	  if(j+1 < im2 && par2[j] != par2[j+1])
	    side2 = 1;
	  else
	    side2 = -1;

	  for(i=0; i<im1; i++)
	    {
	      parvalue[0] = par1[i];
	      if(i+1 < im1 && par1[i] != par1[i+1])
		side1 = 1;
	      else
		side1 = -1;

	      s1425( ps, 1, 1, side1, side2, parvalue, &leftknot1,
			  &leftknot2, derive, &stat );
	      if(stat < 0) goto error;
	      for(k=0; k<idim; k++, index++)
		{
		  ep[index] = derive[k];
		  eder10[index] = derive[k+idim];
		  eder01[index] = derive[k+2*idim];
		  eder11[index] = derive[k+3*idim];
		}
	    }
	}
      if (stat < 0) goto error;

      /* Compute Hermite interpolant */

      s1530(ep, eder10, eder01, eder11, par1, par2, im1, im2,
	    idim, &osurf2, &stat);
      if(stat < 0) goto error;
    }

  /* Perform datareduction on the bicubic surface.
   * Use tolerance error2                         */

  s1965(osurf2, error2, nend, iopen1, iopen2, edgeps, iopt, itmax,
	rs, error1, &stat);
  if(stat < 0) goto error;

  /*
   * Success.
   * --------
   */

  *jstat = 0;
  goto out;


  /*
   * Error in allocation.
   * --------------------
   */

  err101:
    *jstat = -101;
    s6err("s1968",*jstat,kpos);
    goto out;

  /* Error in input */

  err103:
    *jstat = -103;
    s6err("s1968",*jstat,kpos);
    goto out;

  /*
   * Error in input, pointer to SISLSurf was SISL_NULL pointer.
   * -----------------------------------------------------
   */

  err150:
    *jstat = -150;
    s6err("s1968",*jstat,kpos);
    goto out;

  /*
   * Error in lower level routine.
   * -----------------------------
   */

   error:
    *jstat = stat;
    s6err("s1968",*jstat,kpos);
    goto out;

  /*
   * Exit s1968.
   * -----------
   */

  out:
    if(error1 != SISL_NULL)    freearray(error1);
    if(error2 != SISL_NULL)    freearray(error2);
    if(error1472 != SISL_NULL) freearray(error1472);
    if(newcoeff != SISL_NULL)  freearray(newcoeff);
    if(ep != SISL_NULL)        freearray(ep);
    if(eder10 != SISL_NULL)    freearray(eder10);
    if(eder01 != SISL_NULL)    freearray(eder01);
    if(eder11 != SISL_NULL)    freearray(eder11);
    if(derive != SISL_NULL)    freearray(derive);
    if(par1 != SISL_NULL)      freearray(par1);
    if(par2 != SISL_NULL)      freearray(par2);
    if(osurf1 != SISL_NULL)    freeSurf(osurf1);
    if(osurf2 != SISL_NULL)    freeSurf(osurf2);

    return;
}
