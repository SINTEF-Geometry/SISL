/*
 * Copyright (C) 1998, 2000-2007, 2010, 2011, 2012, 2013 SINTEF ICT,
 * Applied Mathematics, Norway.
 *
 * Contact information: E-mail: tor.dokken@sintef.no                      
 * SINTEF ICT, Department of Applied Mathematics,                         
 * P.O. Box 124 Blindern,                                                 
 * 0314 Oslo, Norway.                                                     
 *
 * This file is part of SISL.
 *
 * SISL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Affero General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version. 
 *
 * SISL is distributed in the hope that it will be useful,        
 * but WITHOUT ANY WARRANTY; without even the implied warranty of         
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          
 * GNU Affero General Public License for more details.
 *
 * You should have received a copy of the GNU Affero General Public
 * License along with SISL. If not, see
 * <http://www.gnu.org/licenses/>.
 *
 * In accordance with Section 7(b) of the GNU Affero General Public
 * License, a covered work must retain the producer line in every data
 * file that is created or manipulated using SISL.
 *
 * Other Usage
 * You can be released from the requirements of the license by purchasing
 * a commercial license. Buying such a license is mandatory as soon as you
 * develop commercial activities involving the SISL library without
 * disclosing the source code of your own applications.
 *
 * This file may be used in accordance with the terms contained in a
 * written agreement between you and SINTEF ICT. 
 */

#include "sisl-copyright.h"

/*
 *
 * $Id: sh1834.c,v 1.5 2001-03-19 15:59:05 afr Exp $
 *
 */


#define SH1834

#include "sislP.h"

/*
* Forward declarations.
* ---------------------
*/

#if defined(SISLNEEDPROTOTYPES)
static void sh1834_s9mat2d(double [],double []);
static void sh1834_s9mat3d(double [],double [],double []);
#else
static void sh1834_s9mat2d();
static void sh1834_s9mat3d();
#endif

#if defined(SISLNEEDPROTOTYPES)
void sh1834(SISLObject *po1,SISLObject *po2,double aepsge,int idim,
	    double edir1[],double edir2[],int *jstat)
#else
void sh1834(po1,po2,aepsge,idim,edir1,edir2,jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int    idim;
     double edir1[];
     double edir2[];
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Perform a box-test to check if two objects overlap in
*              the rotated coordinate system where edir1 defines the
*              x-axis and (edir1 x edir2) defines the z-axis if
*              idim = 3.
*
*
*
* INPUT      : po1    - Pointer to the first object.
*              po2    - Pointer to the second object.
*              aepsge - Geometry resolution.
*              idim   - The dimension of the space in which the objects
*                       lie. idim = 2 or idim = 3.
*              edir1  - First direction vector. Defines x-axis in rotated
*                       coordinate system.
*              edir2  - Second direction vector.
*
*
*
* OUTPUT     : jstat  - status messages
*                                 = 2      : Boundaries just touch.
*                                 = 1      : Rotated SISLbox overlaps point.
*                                 = 0      : No overlap.
*                                 < 0      : error
*
*
* METHOD     : The coordinate system is rotated such that if idim = 2,
*              the x-axis of the new coordinate system is parallell to
*              the vector edir1. If idim = 3, the cross-product of edir1
*              and edir2 is rotated to be parallell to the z-axis and
*              edir1 rotated to be parallell to the x-axis. The objects
*              are moved into this rotated coordinate system and a
*              box-test is performed.
*
*
* REFERENCES :
*
*-
* CALLS      : s1790  - Perform box test.
*              s6scpr - Compute scalar product of two vectors.
*              sh1834_s9mat2d - Set up rotation matrix for 2 dimensional space.
*              sh1834_s9mat3d - Set up rotation matrix for 3 dimensional space.
*
* WRITTEN BY : Vibeke Skytt, SI, 88-06.
* REVISED BY : Vibeke Skytt, SI, 91-01.
* REVISED BY : Mike Floater, SI, 94-07. Updated for rationals
*                                        (s1834 didn't need changing).
* REVISED BY : Vibeke Skytt, SINTEF Oslo, 94-11. Introduced 2D point-surface
*                                                intersection.
*
*********************************************************************
*/
{
  int kstat = 0;   /* Local status variable.                     */
  int kpos = 0;    /* Position of error.                         */
  int kinnerexp = 12; /* Expand box in the inner. No rotation.   */
  int kn1=0,kn2=0;     /* Number of coefficients of objects.     */
  double *sc1,*sc2;/* Pointers to coefficients of objects.       */
  double *scoef1=SISL_NULL;  /* Rotated coefficients of first object. */
  double *scoef2=SISL_NULL;  /* Rotated coefficients of second object.*/
  double *smat=SISL_NULL;    /* Rotation matrix.                      */
  double *s1,*s2,*s3,*s4,*s5;  /* Pointers used to traverse arrays. */
  SISLObject *qo1=SISL_NULL; /* First object after rotation.          */
  SISLObject *qo2=SISL_NULL; /* Second object after rotation.         */
  /*  long time_before;
  long time_used=0; */

  double *rc1,*rc2;     /* Pointers to homogeneous coefficients. */
  double *rcoef1=SISL_NULL;  /* Possibly homogeneous coefficients.    */
  double *rcoef2=SISL_NULL;  /* Possibly homogeneous coefficients.    */
  int ikind1=0, ikind2=0;   /* Kinds of objects 1 and 2.         */
  int i,i1,i2,j,k;      /* Loop variables.                       */

  /* Test input.  */

  if (idim != 2 && idim != 3) goto err105;

  /* Fetch coefficients of the objects. */

  if (po1->iobj == SISLCURVE)
  {
     kn1 = po1->c1->in;
     sc1 = po1->c1->ecoef;
     rc1 = po1->c1->rcoef;
     ikind1 = po1->c1->ikind;
  }
  else if (po1->iobj == SISLSURFACE)
  {
     kn1 = po1->s1->in1*po1->s1->in2;
     sc1 = po1->s1->ecoef;
     rc1 = po1->s1->rcoef;
     ikind1 = po1->s1->ikind;
  }
  else
  {
     kn1 = 1;
     sc1 = po1->p1->ecoef;
     rc1 = SISL_NULL;
     ikind1 = 1;
  }

  if (po2->iobj == SISLCURVE)
  {
     kn2 = po2->c1->in;
     sc2 = po2->c1->ecoef;
     rc2 = po2->c1->rcoef;
     ikind2 = po2->c1->ikind;
  }
  else if (po2->iobj == SISLSURFACE)
  {
     kn2 = po2->s1->in1*po2->s1->in2;
     sc2 = po2->s1->ecoef;
     rc2 = po2->s1->rcoef;
     ikind2 = po2->s1->ikind;
  }
  else
  {
     kn2 = 1;
     sc2 = po2->p1->ecoef;
     rc2 = SISL_NULL;
     ikind2 = 1;
  }

  /* Allocate space for local parameters.  */

  if ((scoef1 = newarray(idim*kn1,DOUBLE)) == SISL_NULL) goto err101;
  if ((scoef2 = newarray(idim*kn2,DOUBLE)) == SISL_NULL) goto err101;
  if ((smat = new0array(idim*idim,DOUBLE)) == SISL_NULL) goto err101;

  /* Find the rotation matrix.  */

  if (idim == 2)

    /* After normalization edir1[0] will contain the cosine of the
       rotation angle and edir1[1] will contain the sine.           */

     sh1834_s9mat2d(smat,edir1);
  else

    /* Set up the rotation matrix when idim = 3. (edir1 x edir2) is
       rotated to be parallell to the z-axis and edir1 to be parallell
       to the x-axis.                                                   */

  sh1834_s9mat3d(smat,edir1,edir2);

  /* The objects is moved into the new coordinate system by rotating
     them using the rotation matrix.                                 */

  /* Rotate first object. */

    for (s2=sc1,s4=s2+idim*kn1,s5=scoef1; s2<s4; s2+=idim)
     for (s1=smat,s3=smat+idim*idim; s1<s3; s1+=idim,s5++)
	*s5 = s6scpr(s1,s2,idim);

  /* Rotate second object. */

  for (s2=sc2,s4=s2+idim*kn2,s5=scoef2; s2<s4; s2+=idim)
     for (s1=smat,s3=smat+idim*idim; s1<s3; s1+=idim,s5++)
	*s5 = s6scpr(s1,s2,idim);

  /* Make rotated objects.  */

  if ((qo1 = newObject(po1->iobj)) == SISL_NULL) goto err101;
  if ((qo2 = newObject(po2->iobj)) == SISL_NULL) goto err101;

  if(ikind1 == 2 || ikind1 == 4)
  {
      if ((rcoef1 = newarray((idim+1)*kn1,DOUBLE)) == SISL_NULL) goto err101;
      for(i=0,i1=0,i2=0; i<kn1; i++)
      {
	  k = i1 + idim;
	  for(j=0; j<idim; j++, i1++, i2++)
	  {
	      rcoef1[i1] = scoef1[i2] * rc1[k];
	  }
	  rcoef1[i1] = rc1[k];
	  i1++;
      }
  }
  else
  {
      rcoef1 = scoef1;
  }

  if(ikind2 == 2 || ikind2 == 4)
  {
      if ((rcoef2 = newarray((idim+1)*kn2,DOUBLE)) == SISL_NULL) goto err101;
      for(i=0,i1=0,i2=0; i<kn2; i++)
      {
	  k = i1 + idim;
	  for(j=0; j<idim; j++, i1++, i2++)
	  {
	      rcoef2[i1] = scoef2[i2] * rc2[k];
	  }
	  rcoef2[i1] = rc2[k];
	  i1++;
      }
  }
  else
  {
      rcoef2 = scoef2;
  }


  if (po1->iobj == SISLCURVE)
  {
     if ((qo1->c1 = newCurve(po1->c1->in,po1->c1->ik,po1->c1->et,
			     rcoef1,po1->c1->ikind,idim,0)) == SISL_NULL)
	goto err101;
     /* printf("Rotated box test. Curve - "); */
  }
  else if (po1->iobj == SISLSURFACE)
  {
     if ((qo1->s1 = newSurf(po1->s1->in1,po1->s1->in2,po1->s1->ik1,
			    po1->s1->ik2,po1->s1->et1,po1->s1->et2,
			     rcoef1,po1->s1->ikind,idim,0)) == SISL_NULL)
	goto err101;
     /* printf("Rotated box test. Surface - "); */

  }
  else 
  {
     if ((qo1->p1 = newPoint(rcoef1,idim,0)) == SISL_NULL) goto err101;
  }

  if (po2->iobj == SISLCURVE)
  {
     if ((qo2->c1 = newCurve(po2->c1->in,po2->c1->ik,po2->c1->et,
			     rcoef2,po2->c1->ikind,idim,0)) == SISL_NULL)
	goto err101;
     /* printf("curve. "); */
  }
  else if (po2->iobj == SISLSURFACE)
  {
     if ((qo2->s1 = newSurf(po2->s1->in1,po2->s1->in2,po2->s1->ik1,
			    po2->s1->ik2,po2->s1->et1,po2->s1->et2,
			     rcoef2,po2->s1->ikind,idim,0)) == SISL_NULL)
	goto err101;
     /* printf("surface. "); */
  }
  else 
  {
     if ((qo2->p1 = newPoint(rcoef2,idim,0)) == SISL_NULL) goto err101;
  }

  /* Make box test.  */

  /*  time_before = clock();
  boxrot_nmb++; */
  sh1790(qo1,qo2,kinnerexp,aepsge,&kstat);
  /*  time_used = clock() - time_before;
  boxrot_time += time_used; */
  if (kstat < 0) goto error;
		 /* printf("Status = %d \n",kstat); */

  /* Box-test permformed.  */

  *jstat = kstat;
  goto out;

  /* Error in space allocation.  */

 err101: *jstat = -101;
  s6err("sh1834",*jstat,kpos);
  goto out;

  /* Error in input. Dimension not equal to 2 or 3.  */

 err105: *jstat = -105;
  s6err("sh1834",*jstat,kpos);
  goto out;

  /* Error in lower level routine.  */

  error : *jstat = kstat;
  goto out;

 out:

  /* Free space occupied by local arrays and objects.  */

  if (qo1 != SISL_NULL) freeObject(qo1);
  if (qo2 != SISL_NULL) freeObject(qo2);
  if (rcoef1 != SISL_NULL && rcoef1 != scoef1) freearray(rcoef1);
  if (rcoef2 != SISL_NULL && rcoef2 != scoef2) freearray(rcoef2);
  if (scoef1 != SISL_NULL) freearray(scoef1);
  if (scoef2 != SISL_NULL) freearray(scoef2);
  if (smat != SISL_NULL) free0array(smat);

  return;
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1834_s9mat2d(double emat[],double edir[])
#else
static void sh1834_s9mat2d(emat,edir)
     double emat[];
     double edir[];
#endif
/*
*********************************************************************
*
* PURPOSE    : Set up rotation matrix in two dimensions when the x-axis
*              is supposed to be rotated to be parallell with edir.
*
* INPUT      : edir  - Direction of rotated x-axis.
*
* OUTPUT     : emat  - Rotation matrix. emat is supposed to be
*                      initialized to zero before this routine is entered.
*
*********************************************************************
*/
{
  int kstat = 0;   /* Local status variable.              */
  double tlength;  /* Length of vector edir.              */
  double sdir[2];  /* Normalized vertion of vector edir.  */

  tlength = s6norm(edir,2,sdir,&kstat);
  if (kstat == 0)

    /* Length of edir equal to zero. Let the rotation matrix be
       the identity matrix.                                      */

    emat[0] = emat[3] = (double)1.0;
  else
    {

      /* Make rotation matrix.  */

      emat[0] = sdir[0];
      emat[1] = sdir[1];
      emat[2] = sdir[1];
      emat[3] = -sdir[0];
    }
}

#if defined(SISLNEEDPROTOTYPES)
static void
  sh1834_s9mat3d(double emat[],double edir1[],double edir2[])
#else
static void sh1834_s9mat3d(emat,edir1,edir2)
     double emat[];
     double edir1[];
     double edir2[];
#endif
/*
*********************************************************************
*
* PURPOSE    : Set up rotation matrix in three dimensions when edir1
*              is supposed to be rotated to be parallell with the x-axis
*              and (edir1 x edir2) to be parallell with the z-axis.
*
* INPUT      : edir1 - First direction vector.
*              edir2 - Second direction vector.
*
* OUTPUT     : emat  - Rotation matrix. The matrix is supposed to be
*                      initialized to zero before this routine is enterd.
*
*********************************************************************
*/
{
  int kstat = 0;    /* Local status variable.                         */
  double snorm[3];  /* Cross-product of edir1 and edir2.              */
  double sdir[3];   /* Normalized vertion of edir1.                   */
  double *s1;       /* Pointer into emat array.                       */
  double tleng1,tleng2; /* Length of snorm and edir1 respectively.    */
  double ta1,ta2,ta3,tb1,tb2,tb3,td1,td2,tl1,tl2,tl3; /* Help variables. */

  /* Calculate cross-product of edir1 and edir2.  */

  s6crss(edir1,edir2,snorm);

  /* Normalize snorm.  */

  tleng1 = s6norm(snorm,3,snorm,&kstat);

  /* Normalize edir1.  */

  tleng2 = s6norm(edir1,3,sdir,&kstat);

  /* Initialize help variables.  */

  ta1 = snorm[0];
  ta2 = snorm[1];
  ta3 = snorm[2];
  tl1 = sqrt(ta2*ta2+ta3*ta3);

  /* Set up rotation matrix.  */

  if ((DEQUAL(tleng1,DZERO) || DEQUAL(tl1,DZERO)) && DEQUAL(tleng2,DZERO))

    /* The rotation matrix is the identity matrix.  */

    emat[0] = emat[4] = emat[8] = (double)1.0;
  else if (DEQUAL(tleng1,DZERO) || DEQUAL(tl1,DZERO))
    {

      /* The rotation matrix is supposed to rotate edir1 to be parallell
	 to the x-axis.                                                   */

      tb1 = sdir[0];
      tb2 = sdir[1];
      tb3 = sdir[2];
      tl3 = sqrt(tb1*tb1+tb2*tb2);

      if (DEQUAL(tl3,DZERO)) emat[0] = emat[4] = emat[8] = (double)1.0;
      else
	{
	  s1      = emat;
	  *(s1++) = tb1;
	  *(s1++) = tb2;
	  *(s1++) = tb3;
	  *(s1++) = -tb2/tl3;
	  *(s1++) = tb1/tl3;
	  *(s1++) = DZERO;
	  *(s1++) = -tb1*tb3/tl3;
	  *(s1++) = -tb2*tb3/tl3;
	  *(s1++) = tl3;
	}
    }
  else
    {
      td1 = edir1[0]/tl1;
      td2 = (ta3*edir1[1] - ta2*edir1[2])/tl1;
      tl2 = sqrt(td1*td1+td2*td2);

      if (DEQUAL(tl2,DZERO))
	{

	  /* The normal snorm is rotated to be parallell to the z-axis. */

	  s1      = emat;
	  *(s1++) = tl1;
	  *(s1++) = -ta1*ta2/tl1;
	  *(s1++) = -ta1*ta3/tl1;
	  *(s1++) = DZERO;
	  *(s1++) = ta3/tl1;
	  *(s1++) = -ta2/tl1;
	  *(s1++) = ta1;
	  *(s1++) = ta2;
	  *(s1++) = ta3;
	}
      else
	{

	  /* The normal is rotated to be parallell to the z-axis and edir1
	     to be parallell to the x-axis.                                 */

	  s1      = emat;
	  *(s1++) = td1*tl1/tl2;
	  *(s1++) = (-ta1*ta2*td1 + ta3*td2)/(tl1*tl2);
	  *(s1++) = (-ta1*ta3*td1 - ta2*td2)/(tl1*tl2);
	  *(s1++) = -td2*tl1/tl2;
	  *(s1++) = (ta1*ta2*td2 + ta3*td1)/(tl1*tl2);
	  *(s1++) = (ta1*ta3*td2 - ta2*td1)/(tl1*tl2);
	  *(s1++) = ta1;
	  *(s1++) = ta2;
	  *(s1++) = ta3;
	}
    }
}
