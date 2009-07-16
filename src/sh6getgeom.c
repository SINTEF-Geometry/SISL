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
 * $Id: sh6getgeom.c,v 1.2 2001-03-19 15:59:07 afr Exp $
 *
 */


#define SH6GETGEOM

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6getgeom(SISLObject *ob, int obnr, SISLIntpt *pt,
		 double **geom, double **norm, double aepsge, int *jstat)
#else
void sh6getgeom(ob,obnr,pt,geom,norm,aepsge,jstat)
   SISLObject *ob;
   int        obnr;
   SISLIntpt  *pt;
   double     **geom;
   double     **norm;
   double     aepsge;
   int        *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt, get the geometric data.
*
*
* INPUT      : ob	- Pointer to geometric object.
*	       obnr	- 1 or 2 number of object to evaluate.
*	       pt       - Pointer to the Intpt.
*
*
* OUTPUT     : geom    	- Geometric data for object.
*              geom2    - Normal if existing, otherwise SISL_NULL(dummy).
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. May 91.
* REVICED BY : Arne Laksaa, SI, Oslo, Norway. May 91.
* CORRECTED BY: UJK (Lets have a OK status exit !)
*********************************************************************
*/
{
   int kgeom;	/* Number of doubles pr object describing geometry. */
   int dim;	/* Geometric dimension. */
   int kpar;	/* Index of the parameter value of the object in pt. */
   int kstat;
   int left1=0,left2=0;
   double *val;
   
   /* UJK */
   *jstat = 0;

   kgeom = (obnr == 1 ? pt->size_1 : pt->size_2);
   
   if (ob->iobj == SISLPOINT)      dim = ob->p1->idim;
   else if (ob->iobj == SISLCURVE) dim = ob->c1->idim;
   else if (ob->iobj == SISLSURFACE)  dim = ob->s1->idim;
   
   kpar = (obnr == 1 ? 0 : (pt->ipar - ob->iobj));

   if (!kgeom)
      switch(ob->iobj)
      {
	 case SISLPOINT:
	    (*geom) = ob->p1->ecoef; 
	    (*norm) = SISL_NULL;
            return;	    
	 case SISLCURVE:
	    val = newarray(2*dim,DOUBLE);
	    shevalc(ob->c1,1,pt->epar[kpar],aepsge,&left1,val,&kstat);
	    if (kstat < 0) goto err1;
	    if (obnr == 1)
	    {
	       pt->geo_data_1 = val;
	       pt->size_1 = 2*dim;
	       kgeom = pt->size_1;
	    }
	    else
	    {
	       pt->geo_data_2 = val;
	       pt->size_2 = 2*dim;
	       kgeom = pt->size_2;
	    }
	    
	    break;
	 case SISLSURFACE:
	    val = newarray(7*dim,DOUBLE);
	    s1421(ob->s1,2,pt->epar+kpar,&left1,&left2,val,val+6*dim,&kstat);
	    if (kstat < 0) goto err1;
	    if (obnr == 1)
	    {
	       pt->geo_data_1 = val;
	       pt->size_1 = (dim == 3 ? 7 : 6)*dim;
	       kgeom = pt->size_1;
	    }
	    else
	    {
	       pt->geo_data_2 = val;
	       pt->size_2 = (dim == 3 ? 7 : 6)*dim;
	       kgeom = pt->size_2;
	    }

	    break;
      }
	    
   
   (*geom) = (obnr == 1 ? pt->geo_data_1 : pt->geo_data_2);
   
   if (ob->iobj == SISLSURFACE) (*norm) = (*geom) + kgeom - dim;
   else				(*norm) = SISL_NULL;
   goto out;
   
   err1: *jstat = kstat;
   goto out;
   
   out :
      return;
}
