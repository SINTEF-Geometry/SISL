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
 *
 *
 */


#define S1516

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s1516(double ep[],double epar[],int im,int idim,
           double **ev,int *jstat)
#else
void s1516(ep,epar,im,idim,ev,jstat)
     double ep[];
     double epar[];
     int    im;
     int    idim;
     double **ev;
     int    *jstat;
#endif
/*
************************************************************************
*
* Purpose:   To estimate the first derivative at each point in a sequence.
*
* Input:
*          ep     - Array containing the point in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          epar   - Parametrization array. The array should be increasing
*                   in value.
*          im     - Number of point and derivatives
*          idim   - The dimension of the space the points and derivatives
*                   lie in
* Output:
*          ev     - Pointer to array containing the derivatives in sequence
*                   (x,y,..,x,y,..), length idim * im.
*          jstat  - Status variable
*                    < 0 - Error.
* Method:
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Michael Floater, SI 1993-10
*
*********************************************************************
*/
{
  int ki,kj;          /* Loop variables                              */
  int kk;             /* Polynomial order                            */
  int kpos=0;         /* Position of error                           */
  int kstat=0;        /* Status variable                             */
  double *gpar;
  int kcnsta;
  int kcnend;
  int iopen;
  int iorder;
  int ileft;
  double *ntype;
  SISLCurve *qc;
  int knbpar;
  double *evtemp;
  double cendpar;
  double *eder;




  /* Check input */

  if (idim < 1 || im < 2) goto err102;


  /* Allocate array for derivatives */

  evtemp    = newarray(idim*im,DOUBLE);
  if (evtemp == SISL_NULL) goto err101;

  ntype    = newarray(im,DOUBLE);
  if (ntype == SISL_NULL) goto err101;

  for(ki=0; ki<im; ki++)
  {
      ntype[ki] = 1.0;
  }

  eder    = newarray(2 * idim,DOUBLE);
  if (eder == SISL_NULL) goto err101;



  kcnsta = 1;
  kcnend = 1;
  iopen = 1;
  iorder = 4;

  s1358(ep, im, idim, ntype, epar, kcnsta, kcnend, iopen, iorder,
        epar[0],&cendpar, &qc, &gpar, &knbpar, &kstat);
     if(kstat < 0) goto error;

  for(ki=0,kk=0; ki<im; ki++,kk+=idim)
  {
      s1221(qc,1,epar[ki],&ileft,eder,&kstat);
      if(kstat < 0) goto error;

      for(kj=0; kj<idim; kj++)
      {
          evtemp[kk+kj] = eder[idim+kj];
      }
  }


  /* Calculation completed */

  /* Set result. */

  (*ev) = evtemp;

  *jstat = 0;
  goto out;



  /* Error in space allocation */

 err101: *jstat = -101;
  s6err("s1516",*jstat,kpos);
  goto out;


  /* Error in input. */

 err102: *jstat = -102;
  s6err("s1516",*jstat,kpos);
  goto out;

  /* Error in lower level routine. */

 error:  *jstat =kstat;
  s6err("s1516",*jstat,kpos);
  goto out;

 out:
  if (ntype != SISL_NULL) freearray(ntype);
  if (gpar != SISL_NULL) freearray(gpar);
  if (eder != SISL_NULL) freearray(eder);

  return;
}
