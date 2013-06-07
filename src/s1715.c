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
 * $Id: s1715.c,v 1.2 2001-03-19 15:58:52 afr Exp $
 *
 */


#define S1715

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1715(SISLCurve *pc1,SISLCurve *pc2,int iend1,int iend2,SISLCurve **rcnew,int *jstat)
#else
void s1715(pc1,pc2,iend1,iend2,rcnew,jstat)
     SISLCurve *pc1;
     SISLCurve *pc2;
     int   iend1;
     int   iend2;
     SISLCurve **rcnew;
     int   *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : To join one end of one B-spline curve with one end of
*              another B-spline curve by translating the second curve.
*              If pc1 is to be joined at the start the direction of the
*              curve is turned, and if pc2 is to be joined at the end
*              the direction of this curve is turned. This means that
*              pc1 always is at the beginning at the new curve.
*
*
*
* INPUT      : pc1     - First curve to join.
*              pc2     - Second curv to join.
*              iend1   - True(1) if the first curv is to be
*                        joined et the end, else false(0).
*              iend2   - True(1) if the second curv is to be
*                        joined at the end, else false(0).
*
*
*
* OUTPUT     : rcnew   - The new joined curve.
*              jstat   - status messages
*                                         > 0      : warning
*                                         = 0      : ok
*                                         < 0      : error
*
*
* METHOD     : We have to turn the first curve if iend1=0
*              and the second curve if iend2=1.
*              Then we translate the second curve so the gap
*              between the first and second curvde disappear.
*
*
* REFERENCES :
*
*-
* CALLS      : newCurve  - Allocate space for a new curve-object.
*              freeCurve - Free space occupied by given curve-object.
*              S1750.C   - Making one curve at a higer degree.
*
* WRITTEN BY : Arne Laksaa, SI, 88-06.
* REVISED BY : Johannes Kaasa, SI, Sep 1991 (Introduced NURBS)
*
**********************************************************************/
{
  int kstat=0;            /* Local status variable.                     */
  int kpos=0;             /* Position of error.                         */
  int kcopy=0;            /* To mark if pc1 (1) or pc2 (2) is
			     changed to point at a local copy.          */
  int km1=0,km2=0;        /* Knot mutiplicety at the end to join.       */
  int km2end=0;           /* Knot mutiplicety at the end of the
			     second curve.                              */
  int kk;                 /* Order of the new curve.                    */
  int kn;                 /* Number of the vertices in the new curve.   */
  int kdim;               /* Dimensjon of the space in whice curve lies.*/
  int routdim;            /* Rational dimension of the output curve.    */
  int kn1=pc1->in;        /* Number of vertices in the old curves.      */
  int kn2=pc2->in;        /* Number of vertices in the old curves.      */
  int ki,kj;              /* Control variable in loop, and others.      */
  double tdel;            /* The translation of the knots to the
			     second curve.                              */
  double *s1,*s2,*s3; 	  /* Pointers used in loop.                     */
  double *stran=SISL_NULL;     /* The translation vector to vertices.        */
  double *st=SISL_NULL;        /* The new knot-vector.                       */
  double *scoef=SISL_NULL;     /* The new vertice.                           */
  SISLCurve *qc=SISL_NULL;     /* Pointer to the new curve-object.           */
  
  int ktype;              /* Type of curves:                            */
                          /* = 1 : Both are B-splines                   */
                          /* = 2 : pc1 is B-spline and pc2 is NURBS     */
                          /* = 3 : pc1 is NURBS and pc2 is B-spline     */
                          /* = 4 : Both are NURBS                       */
  int knumb;              /* Number of vertices.                        */
  double weight;          /* Rational weight.                           */
  double *u1,*u2;         /* Utility pointers into the vertices.        */

  /* Check that we have curves to join. */
  
  if (!pc1 || !pc2) goto err150;

  /* Check that The curves is in the same room, have the same kdim. */
  
  if (pc1->idim != pc2->idim) goto err106;
  else kdim = pc1->idim;

  /* Check the type of the curves. */

  if (pc1->ikind == 2 || pc1->ikind == 4)
    {
      if (pc2->ikind == 2 || pc2->ikind == 4)
        {
          ktype = 4;
          routdim = kdim + 1;
        }
      else
        {
          ktype = 3;
          routdim = kdim + 1;
        }
    }
  else
    {
      if (pc2->ikind == 2 || pc2->ikind == 4)
        {
          ktype = 2;
          routdim = kdim + 1;
        }
      else
        {
          ktype = 1;
          routdim = kdim;
        }
    }
  
  /* Allocate a kdim array to store the translation of the second curv.*/
  
  if ((stran=newarray(kdim,double)) == SISL_NULL) goto err153;
  
  /* Checking the order of the curves, and raise the order if nessesary.*/
  
  if (pc1->ik < pc2->ik)
    {
      kcopy=1;
      kk=pc2->ik;
      s1750(pc1,kk,&pc1,&kstat);
      if (kstat) goto err153;
    } 
  else
    if (pc2->ik < pc1->ik)
      {
	kcopy=2;
	kk=pc1->ik;
	s1750(pc2,kk,&pc2,&kstat);
	if (kstat) goto err153;
      } 
    else
      kk = pc1->ik;
  
  /* Finding the knot multiplicity at the juinction, km1 km2.
     At the end thats  going to be the end of the new curve
     we also need to know the knot mutiticeply, km2end. */
  
  /* Having raised the order of the curves if necessary,
     remember the number of vertices in the two curves. */

  kn1=pc1->in;
  kn2=pc2->in;

  if (iend1) 
    while (pc1->et[kn1+kk-1-km1] == pc1->et[kn1+kk-1]) km1++;
  else        
    while (pc1->et[km1] == *pc1->et) km1++;
  if (iend2)
    {
      while (pc2->et[kn2+kk-1-km2] == pc2->et[kn2+kk-1]) km2++;
      while (pc2->et[km2end] == *pc2->et) km2end++;
    } 
  else
    {
      while (pc2->et[km2] == *pc2->et) km2++;
      while (pc2->et[kn2+kk-1-km2end] == pc2->et[kn2+kk-1]) km2end++;
    }
  
  /* Find the number of vertices in the new curve. */
  
  kn = kn1 + kn2 + 3*kk - km1 - km2 - km2end -1;
  
  /* Allocating the new arrays to the new curve. */
  
  if ((st=newarray(kn+kk,double))==SISL_NULL) goto err101;
  if ((scoef=newarray(kn*routdim,double))==SISL_NULL) goto err101;
  
  /* Copying the knotvectors from the old curve to the new curves */
  /****************************************************************/
  
  /* The first curve. */
  
  if (iend1)     /* The junction is at the end of the first curve. */
    {
      /* Copying all knots from the first curve that is different
	 from the last knot. */
      
      memcopy(st,pc1->et,kn1+kk-km1,double);
      
      /* Making a kk-1 touple knot at the junction. */
      
      for (s1=st+kn1+kk-km1,s2=s1+kk-1; s1<s2; s1++)
	*s1=pc1->et[kn1+kk-1];
      
    } 
  else     /* The junction is at the beginning of the first curve. */
    {
      /* Computing the factor to turn the first knotvector. */
      
      tdel = *pc1->et + pc1->et[kn1+kk-1];
	
      /* Copying and turning the first knot vector except the
	 knots at the junction. */
      
      for (s1=st,s2=pc1->et+km1,s3=pc1->et+kn1+kk-1; s2<=s3; s1++,s3--)
	*s1 = tdel - *s3;
      
      /* Making a kk-1 touple knot at the junction. */
      
      for (s2--,s3=s1+kk-1; s1<s3; s1++) *s1= tdel - *s2;
    }
  
  /* The second curve. */
  
  s2 = st+kn+kk-max(0,kk-km2end); /* The border for the last exsisting knot 
				     in the new knot vector. */
  
  if (!iend2)  /* The junction is at the begining of the second curve. */
    {
      /* Computing what the second knot vector has to be translated
	 to get a kontinue total knotvector. */
      
      tdel = s1[-1] - *pc2->et;
      
      /* copying and translating all knots except the knots
	 at the junction. */
      
      for (s3=pc2->et+km2; s1<s2; s1++,s3++) *s1 = *s3 + tdel;
    } 
  else
    {
      /* Coputing a factor to both translate and turn the knots. */
      
      tdel = pc2->et[kn2+kk-1] + s1[-1];
      
      /* Turning and translating all knots exept the knots
	 at the junction. */
	
      for (s3=pc2->et+kn2+kk-km2-1; s1<s2; s1++,s3--) *s1 =tdel- *s3;
    }
  
  /* Inserting new knots such that we have a kk touple knot at the end.*/
  
  for (ki=0; ki<kk-km2end; ki++)  s1[ki] = s1[-1];
  
  /* Copying the coeffesientvectors to the new curves.*/
  /***************************************************/
  
  /* Copying the first coeffisientvector. */
  
  knumb = min(kn1,kn1+kk-km1);
  ki = routdim*knumb;
  if (iend1)                        /* Just copying. */
    {
      if (ktype == 1)
        memcopy(scoef,pc1->ecoef,ki,double);
      else if (ktype == 2)
        {
          for (kj=0; kj<knumb; kj++)
            {
              u1 = scoef + kj*routdim;
              u2 = pc1->ecoef + kj*kdim;
              memcopy(u1,u2,kdim,double);
              scoef[kj*routdim + kdim] = 1.;
            }
        }
      else if (ktype == 3 || ktype == 4) 
        memcopy(scoef,pc1->rcoef,ki,double);
      s1 = scoef +ki;
    }
  else                              /* Copying back to front. */
    {
      if (ktype == 2)
        {
          s2 = pc1->ecoef;
          s3 = s2 + kdim*(knumb - 1);
          for (s1=scoef; s2<=s3; s3-=2*kdim)
            {
              for (ki=0; ki<kdim; ki++,s1++,s3++)  *s1 = *s3;
              *s1 = 1.;
              s1++;
            } 
        }
      else
        {
          if (ktype == 1)
            s2 = pc1->ecoef;
          else
            s2 = pc1->rcoef; 
          for (s1=scoef,s3=s2+ki-routdim; s2<=s3; s3-=2*routdim)
            for (ki=0; ki<routdim; ki++,s1++,s3++)  *s1 = *s3;
        }
    }

  /* If there is less than a kk touple knot at the end of the first curve
     than we have to inserte zeroes. */
  
  for (s2=s1+routdim*max(0,kk-km1); s1<s2; s1+=routdim)
    {
      for (ki=0; ki<kdim; ki++)
        s1[ki] = DZERO;
      if (ktype != 1)
        s1[kdim] = 1.;
    }
  
  /* Compute the translation of the second curv. */
  
  for (ki=0; ki<kdim; ki++)
    {
      if (km2<kk) stran[ki] = DZERO;
      else
	stran[ki] = iend2? pc2->ecoef[kdim*(kn2-max(0,km2-kk)-1)+ki]:
	  pc2->ecoef[kdim*max(0,km2-kk)+ki];
      if (km1>=kk)
	stran[ki] -= iend1? pc1->ecoef[kdim*(kn1-max(0,km1-kk)-1)+ki]:
	  pc1->ecoef[kdim*max(0,km1-kk)+ki];
    }
  
  
  /* Copying the second coefficientvector. */
  
  /* Findig the startpoint for copying in the old coeffisient vector. */
  
  if (iend2)
    {
      if (ktype == 1 || ktype == 3)
        s3=pc2->ecoef+kdim*(kn2-max(0,km2-kk)-1);
      else
        s3=pc2->rcoef+routdim*(kn2-max(0,km2-kk)-1);
    }
  else
    {
      if (ktype == 1 || ktype == 3) 
        s3=pc2->ecoef+kdim*max(0,km2-kk);
      else
        s3=pc2->rcoef+routdim*max(0,km2-kk);
    }
  
  if (km2<kk)
    {
      /* If km2<kk-1 we have to insert zeroes and than transform. */
      
      for (kj=km2+1; kj<kk; kj++,s2+=routdim)
        {
	  for (ki=0; ki<kdim; ki++,s1++) *s1 = -stran[ki];
          if (ktype != 1)
            {
              *s1 = 1.;
              s1++;
            }
        }       

      /* Copying and transforming the first coeffisients. */
      
      if (ktype == 1) 
        for (ki=0; ki<kdim; ki++,s1++,s3++) *s1 = *s3 - stran[ki];
      else
        {
          if (ktype == 3)
            weight = 1.;
          else
            weight = s3[kdim];
          for (ki=0; ki<kdim; ki++,s1++,s3++) *s1 = *s3 - stran[ki]*weight;
          *s1 = weight;
          s1++;
          if (ktype != 3) s3++;
        } 
    } 
  else
    if (ktype == 1 || ktype == 3)
      s3+=kdim;     /* Skiping the first coeffisients. */
    else  
      s3+=routdim;     /* Skiping the first coeffisients. */

  /* Copying and transforming the coeffisient from the second curve. */
  
  for (s2=scoef+routdim*min(kn,kn-kk+km2end); s1<s2;)
    {
      if (iend2)
        {
          if (ktype == 1 || ktype == 3)
            s3-=2*kdim;
          else
            s3-=2*routdim;
        } 
      if (ktype == 1) 
        for (ki=0; ki<kdim; ki++,s1++,s3++) *s1 = *s3 - stran[ki];
      else
        {
          if (ktype == 3)
            weight = 1.;
          else
            weight = s3[kdim];
          for (ki=0; ki<kdim; ki++,s1++,s3++) *s1 = *s3 - stran[ki]*weight;
          *s1 = weight;
          s1++;
          if (ktype != 3) s3++;
        }
    }
  
  /* Insert and transform from zeroes if we do not have a kk touple
     knot at the end of the second curve. */
  
  for (kj=0; kj<kk-km2end; kj++)
    {
      for (ki=0; ki<kdim; ki++,s1++)  *s1 = -stran[ki];
      if (ktype != 1)
        {
          *s1 = 1.;
          s1++;
        }
    }
  
  /* Create the new curve. */
  
  if (ktype == 1)
    {
      if ((qc=newCurve(kn,kk,st,scoef,1,kdim,2)) == SISL_NULL) goto err101;
    }
  else
      if ((qc=newCurve(kn,kk,st,scoef,2,kdim,2)) == SISL_NULL) goto err101;
  
  /* Updating output. */
  
  *rcnew = qc;
  *jstat = 0;
  goto out;
  
  
  /* Error. Subrutine error. */
  
 err153:
  *jstat = kstat;
  goto outfree;
  
  
  /* Error. No curve to subdevice.  */
  
 err150:
  *jstat = -150;
  s6err("s1715",*jstat,kpos);
  goto out;
  
  
  /* Error. Different dimensjon of the room.  */
  
 err106:
  *jstat = -106;
  s6err("s1715",*jstat,kpos);
  goto out;
  
  
  /* Error. Allocation error, not enough memory.  */
  
 err101:
  *jstat = -101;
  s6err("s1715",*jstat,kpos);
  goto outfree;
  
  
 outfree:
  if(qc) 
    freeCurve(qc);
  else
    {
      if (st) freearray(st);
      if (scoef) freearray(scoef);
    }
  
  /* Free local used memory. */
  
 out: 
  if (stran) 
    freearray(stran);
  if (kcopy == 1) 
    freeCurve(pc1);
  else
    if (kcopy == 2) freeCurve(pc2);
  return;
}

