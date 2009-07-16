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
 * $Id: s1190.c,v 1.2 2001-03-19 15:58:42 afr Exp $
 *
 */


#define S1190

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1190(SISLObject *po1, double *cmax, double aepsge,int *jstat)
#else
void s1190(po1,cmax,aepsge,jstat)
     SISLObject *po1;
     double        *cmax;
     double        aepsge;
     int           *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Perform a box-test on the one-dimensional object po1 
*              versus the level value cmax;
*
*
* INPUT      : po1    - The object.
*              cmax   - The level value.
*	       aepsge - Geometry resolution.
*                                                                     
*
* OUTPUT     : jstat  - status messages  

*                            = 3      : Only touching in corner(s).
*                            = 2      : The object is constant.
*                            = 1      : No overlap.
*			     = 0      : Boxmax  > level value.
*                            < 0      : error
*
*
* METHOD     :
*
*
* REFERENCES :
*
*-
* CALLS      :
*
* WRITTEN BY : Ulf J. Krystad, SI, 89-05.
*
*********************************************************************
*/                                     
{
  int kstat = 0;        /* Local status error.                        */
  int kpos = 0;         /* Position of error.                         */
  int kcorn = 0;        /* Number of corners in object.               */
  int li[4];	        /* Contains the indexes of the corners.       */
  int kbez = 0;         /* Flag to mark bezier curve or patch.        */
  int kdim;	        /* Dimension of space.			      */
  int in1,in2;	        /* Local number of vertices.     	      */
  int i1;	        /* Counter.                     	      */
  int kmax;             /* Index for the largest value of 
			   the vertices of the object*/
  double *tmin1,*tmax1; /* Smallest and largest value of 
			   the vertices of the object*/
  double scorn[4];      /* The corner values of the object*/
  
  
  *jstat = 0;  
  
  /* Check kind of first object. */
  
  if (po1->iobj == SISLPOINT)
    {
      kcorn = 0;
      /* Fetch dimention of the object. */
      
      if((kdim = po1->p1->idim) != 1) goto err105;;
      
      
      /* Check if the SISLbox have been computed. */
      
      if (po1->p1->pbox == SISL_NULL)
	{
	  /* If not compute a box. */
	  
	  s1192(po1,aepsge,&kstat);
	  if (kstat<0) goto error;
	}
      
      /* Fetch the SISLbox boarder. */
      
      kmax  = po1->p1->pbox->imax;      
      tmax1 = po1->p1->pbox->emax;
      tmin1 = po1->p1->pbox->emin;
    }
  else
    if (po1->iobj == SISLCURVE)
      {
	/* Fetch dimention of the object. */
	
	if((kdim = po1->c1->idim) != 1) goto err105;;
	
	/* Fetch corners. */
	
	kcorn = 2;
	li[0] = 0;
	li[1] = po1->c1->in - 1;
	scorn[0] = po1->c1->ecoef[li[0]];
	scorn[1] = po1->c1->ecoef[li[1]];
	
	/* Check if we have a bezier curve. */
	
	if (po1->c1->in == po1->c1->ik) kbez = 1;
	
	/* Check if the SISLbox have been computed. */
	
	if (po1->c1->pbox == SISL_NULL)
	  {
	    /* If not compute a box. */
	    
	    s1192(po1,aepsge,&kstat);
	    if (kstat<0) goto error;
	  }
	
	/* Fetch the SISLbox boarder. */
	kmax  = po1->c1->pbox->imax;      	
	tmax1 = po1->c1->pbox->emax;
	tmin1 = po1->c1->pbox->emin;
      }
    else
      if (po1->iobj == SISLSURFACE)
	{
	  /* Fetch dimention of the object. */
	  
	  if((kdim = po1->s1->idim) != 1) goto err105;;
	  
	  
	  kcorn = 4;
	  in1   = po1->s1->in1;
	  in2   = po1->s1->in2;
	  li[0] = 0;
	  li[1] = in1 - 1;
	  li[2] = in1*(in2 - 1);
	  li[3] = in1*in2-1;
	  scorn[0] = po1->s1->ecoef[li[0]];
	  scorn[1] = po1->s1->ecoef[li[1]];
	  scorn[2] = po1->s1->ecoef[li[2]];
	  scorn[3] = po1->s1->ecoef[li[3]];
	  
	  /* Check if we have a bezier patch. */
	  
	  if (po1->s1->in1 == po1->s1->ik1 &&
	      po1->s1->in2 == po1->s1->ik2)    kbez = 1;
	  
	  /* Check if the SISLbox have been computed. */
	  
	  if (po1->s1->pbox == SISL_NULL)
	    {
	      /* If not compute a box. */
	      
	      s1192(po1,aepsge,&kstat);
	      if (kstat<0) goto error;
	    }
	  
	  /* Fetch the SISLbox boarder. */
	  kmax  = po1->s1->pbox->imax;	  
	  tmax1 = po1->s1->pbox->emax;
	  tmin1 = po1->s1->pbox->emin;
	}
      else  goto err121;
  
  
  /* Now we've got the box, do the test: */
  
  if (*cmax - *tmax1 > aepsge)
    *jstat = 1;         /* The object is beyond level value. */
  
  else if (*tmax1 - *tmin1 < aepsge)
    *jstat = 2;         /* The object is of constant value. */ 
  
  else
    /* if (kbez)*/
    {
      /*check for corner max. */
      for (i1=0; i1<kcorn; i1++)
	if (fabs(scorn[i1] - *tmax1) < aepsge)
	  {
	    
	    *jstat = 3;         /* Only corner touching possible.*/
	    break;
	    
	  }
    }
  /*  
    else */
  /*check for absolute corner max. */
  /*    for (i1=0; i1<kcorn; i1++)
	if (kmax == li[i1])
	{
	
	*jstat = 3;    */     /* Only corner touching possible.*/
  /*	  break;
	  
	  }
	  
	  */  
  
  /* Box-test performed. */
  goto out;

  /* Dimensions not equal one. */
  
 err105: *jstat = -105;
  s6err("s1190",*jstat,kpos);
  goto out;
  
  /* Kind of object does not exist. */
  
 err121: *jstat = -121;
  s6err("s1190",*jstat,kpos);
  goto out;
  
  /* Error in lower level routine. */
  
 error:  *jstat = kstat;
  s6err("s1190",*jstat,kpos);
  goto out;
  
 out:	return;
}





