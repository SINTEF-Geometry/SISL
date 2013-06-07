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
 * $Id: s1741.c,v 1.3 2001-03-19 15:58:53 afr Exp $
 *
 */


#define S1741

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void
     s1741(SISLObject *po1,SISLObject *po2,double aepsge,int *jstat)
#else
void s1741(po1,po2,aepsge,jstat)
     SISLObject *po1;
     SISLObject *po2;
     double aepsge;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Check if a object-object intersection is a simple case,
*              i.e. the intersection will result in one single point.
*
*
*
* INPUT      : po1    - First curve in the intersection problem.
*              po2    - Second curve in the intersection problem.
*              aepsge - Geometry resolution.
*
*
* OUTPUT     : jstat  - status messages  
*                                         = 1      : simpel case.
*                                         = 0      : not simpel case.
*                                         < 0      : error.
*
*
* METHOD     :
*
*
* REFERENCES :
*
* CALLS      : s1990   - Making the direction cone of surface.
*              s1991   - Making the direction cone of curve.
*              sh1993  - Simple case of curve in one dimention.
*              sh1994  - Simple case of surface in one dimention.
*              s6err   - Gives error message.
*
* WRITTEN BY : Arne Laksaa, SI, 89-05.
* REVISED BY : Michael Floater, SI, May 1997 to agree with the HP version.
*
*********************************************************************
*/
{
  int kstat = 0;    /* Local status variable.                          */
  int kpos = 0;     /* Position of the error.                          */
  int k1;           /* Control variable in loop.		       */
  double tang;	    /* Angel between two vectors.		       */
  double small_tang;/* Smallest angle between two vectors.	       */
  
  if (po1->iobj == SISLPOINT || po2->iobj == SISLPOINT)
    {
      SISLObject *qo1,*qo2;
      
      if(po1->iobj == SISLPOINT)
	{
	  qo1 = po1;
	  qo2 = po2;
	}
      else
	{
	  qo1 = po2;
	  qo2 = po1;
	}
      
      if (qo2->iobj == SISLCURVE)
	{
	  /* Test if the curve lies in the same space as the point.  */
	  
	  if (qo1->p1->idim != qo2->c1->idim) goto err106;
	  
	  if (qo2->c1->idim == 1)
	    {
	      sh1993(qo2->c1,aepsge,&kstat);
	      
	      *jstat = kstat;
	      goto out;
	    }
	  
	  /* Computing the direction cone of the curve. If the curve
	     have cones greater then pi we just return not a simple case.  */
	  
	  s1991(qo2->c1,aepsge,&kstat);
	  if (kstat < 0) goto error;
	  else if (qo2->c1->pdir->igtpi != 0) goto out2;/* Not a simple case.*/
	  
	  
	  /* Performing a simple case check. */
	  
	  if (qo2->c1->pdir->aang<PIHALF)
	    {
	      /* A simpel case. The iteration is able to
		 find intersection.*/
	      
	      *jstat = 1;
	      goto out;
	    }
	}
      else if (qo2->iobj == SISLSURFACE)
	{
	  /* Test if the surface lies in the same space as the point.  */
	  
	  if (qo1->p1->idim != qo2->s1->idim) goto err106;
	  
	  
	  if (qo2->s1->idim == 1)
	    {
	      sh1994(qo2->s1,aepsge,&kstat);
	      
	      *jstat = kstat;
	      goto out;
	    }
	  else
	    {
	      /* Computing the direction cone of the surface. If the surface
		 have cones greater then pi we just return not a simple case.*/
	      
	      s1990(qo2->s1,aepsge,&kstat);
	      if (kstat < 0) goto error;
	      else if (qo2->s1->pdir->igtpi != 0) goto out2; /*No simple case*/
	      
	      
	      /* Performing a simple case check. */
	      
	      if (qo2->s1->pdir->aang<PIHALF)
		{
		  /* A simpel case. The iteration is able to
		     find intersection.*/
		  
		  
		  *jstat = 1;
		  goto out;
		}
	    }
	}
    }
  else if (po1->iobj == SISLCURVE && po2->iobj == SISLCURVE)
    {
      /* Test if the curves lies in the same space.  */
      
      if (po2->c1->idim != po1->c1->idim) goto err106;
      
      
      
      /* Computing the direction cone of the two curves. If one of them
	 have cones greater then pi we just return not a simple case.  */
      
      s1991(po1->c1,aepsge,&kstat);
      if (kstat < 0) goto error;

      s1991(po2->c1,aepsge,&kstat);
      if (kstat < 0) goto error;

      if (po1->c1->pdir->igtpi != 0) goto out2;  /* Not a simple case.*/
      if (po2->c1->pdir->igtpi != 0) goto out2;  /* Not a simple case.*/
      
      
      /* Computing the angle beetween the senters of the two cones. */
      
      for (tang=DZERO,k1=0;k1<po1->c1->idim;k1++)
	tang += po1->c1->pdir->ecoef[k1]*po2->c1->pdir->ecoef[k1];
      
      if (tang >= DZERO)  tang = min((double)1.0,tang);
      else                tang = max((double)-1.0,tang);
      
      tang = acos(tang);
      
      if (tang > PIHALF)
         small_tang = PI - tang;
      else
         small_tang = tang;
      
      /* Performing a simple case check. */
      
      if ((tang+po1->c1->pdir->aang+po2->c1->pdir->aang)<PI &&
	  (po1->c1->pdir->aang+po2->c1->pdir->aang)<tang)
	{
	  /* A simpel case. The two cones and their mirrors
	     are not intersecting.*/
	  
	  *jstat = 1;
	  goto out;
	}
      else if (po1->c1->idim == 2)
        {
	  *jstat = 0;
	  goto out;
	}
      else if (tang < PI - ANGULAR_TOLERANCE && 
	       tang > ANGULAR_TOLERANCE      &&
	       po1->c1->pdir->aang <= (double)1.3*small_tang &&
	       po2->c1->pdir->aang <= (double)1.3*small_tang)
	 /*po1->c1->pdir->aang <= (double)1.3*tang &&
	       po2->c1->pdir->aang <= (double)1.3*tang)*/
	{
	  s1796(po1->c1,po2->c1,aepsge,tang,&kstat);
	  if (kstat<0) goto error;
	  else *jstat = kstat;
	  goto out;
	}
    }
  else if (po1->iobj == SISLSURFACE && po2->iobj == SISLSURFACE)
    {
      
      /* Test if the surfaces lies in the same space.  */
      
      if (po2->s1->idim != po1->s1->idim) goto err106;
      
      
      
      /* Computing the direction cone of the two surfaces. If one of them
	 have cones greater then pi we just return not a simple case.  */
      
      s1990(po1->s1,aepsge,&kstat);
      if (kstat < 0) goto error;
      
      s1990(po2->s1,aepsge,&kstat);
      if (kstat < 0) goto error;

      if (po1->s1->pdir->igtpi != 0) goto out2;  /* Not a simple case.  */

      if (po2->s1->pdir->igtpi != 0) goto out2;  /* Not a simple case.  */
      
      /* Computing the angle beetween the senters of the two cones. */
      
      for (tang=DZERO,k1=0;k1<po1->s1->idim;k1++)
	tang += po1->s1->pdir->ecoef[k1]*po2->s1->pdir->ecoef[k1];
      
      if (tang >= DZERO)  tang = min((double)1.0,tang);
      else                tang = max((double)-1.0,tang);
      
      tang = acos(tang);
      
      
      /* Performing a simple case check. */
      
      if ((tang+po1->s1->pdir->aang+po2->s1->pdir->aang)<PI &&
	  (po1->s1->pdir->aang+po2->s1->pdir->aang)<tang)
	{
	  /* A simpel case. The two cones and their mirrors
	     are not intersecting.*/
	  
	  po1->psimple = po2;
	  *jstat = 1;
	  goto out;
	}
      else if (tang < PI - ANGULAR_TOLERANCE && 
	       tang > ANGULAR_TOLERANCE      &&
	       po1->s1->pdir->aang <= (double)1.3*tang &&
	       po2->s1->pdir->aang <= (double)1.3*tang)
	{
	  s1795(po1->s1,po2->s1,aepsge,tang,&kstat);
	  if (kstat < 0) goto error;
	  if (kstat == 1) po1->psimple = po2;
	  *jstat = kstat;
	  goto out;
	}
    }
  else if (po1->iobj == SISLCURVE || po2->iobj == SISLCURVE)
    {
      SISLObject *qo1,*qo2;
      
      if(po1->iobj == SISLCURVE)
	{
	  qo1 = po1;
	  qo2 = po2;
	}
      else
	{
	  qo1 = po2;
	  qo2 = po1;
	}
      
      
      /* Test if the surface and curve lies in the same space.  */
      
      if (qo2->s1->idim != qo1->c1->idim) goto err106;
      
      
      
      /* Computing the direction cone of the curve and the surface. If one of
	 them have cones greater then pi we just return not a simple case. */
      
      
      s1990(qo2->s1,aepsge,&kstat);
      if (kstat < 0) goto error;
      
      s1991(qo1->c1,aepsge,&kstat);
      if (kstat < 0) goto error;

      if (qo1->c1->pdir->igtpi != 0) goto out2;  /* Not a simple case.  */
      if (qo2->s1->pdir->igtpi != 0) goto out2;  /* Not a simple case.  */

      
      
      /* Computing the angle beetween the senters of the two cones. */
      
      for (tang=DZERO,k1=0;k1<qo2->s1->idim;k1++)
	tang += qo2->s1->pdir->ecoef[k1]*qo1->c1->pdir->ecoef[k1];
      
      if (tang >= DZERO) tang = min((double)1.0,tang);
      else               tang = max((double)-1.0,tang);
      
      tang = acos(tang);
      
      
      /* Performing a simple case check. */
      
      if (((tang + qo1->c1->pdir->aang) < (PIHALF - qo2->s1->pdir->aang)) ||
	  ((tang - PIHALF - qo1->c1->pdir->aang) > qo2->s1->pdir->aang)) 
	{
	  /* A simpel case. The curve cone or the mirror cone
	     are tottally inside the inverted surface cone. */
	  
	  *jstat = 1;
	  goto out;
	}
      else if (tang < PI - ANGULAR_TOLERANCE && 
	       tang > ANGULAR_TOLERANCE      &&
	       min(tang,fabs(PI-tang)) < 
	       (double)0.8*(PIHALF - qo2->s1->pdir->aang) &&
	       qo1->c1->pdir->aang < (double)0.8*(PIHALF-qo2->s1->pdir->aang))
	{
	  s1797(qo2->s1,qo1->c1,aepsge,tang,&kstat);
	  if (kstat<0) goto error;
	  else *jstat = kstat;
	  goto out;
	}
    }
  

/* Not a simple case. */

out2:	*jstat = 0;
	goto out;

/* Error. Dimensions conflicting.  */

err106: *jstat = -106;
        s6err("s1741",*jstat,kpos);
        goto out;

/* Error in lower level routine.  */

error : *jstat = kstat;
        s6err("s1741",*jstat,kpos);
        goto out;

out:  ;
}
