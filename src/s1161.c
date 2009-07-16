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
 * $Id: s1161.c,v 1.2 2001-03-19 15:58:41 afr Exp $
 *
 */


#define S1161

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s1161(SISLObject *po1,double *cmax,double aepsge,SISLIntdat **pintdat,int *jstat)
#else
void s1161(po1,cmax,aepsge,pintdat,jstat)
     SISLObject *po1;
     double *cmax;
     double aepsge;
     SISLIntdat **pintdat;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : Find all maximum points of an onedimentional object (of type 
*              point, curve or surface). The maximum points found has to
*              be greater or equal to the level value cmax . 
*              In this rouine the outer edges/endpoints of the objects
*               are treated.
*
*
*
* INPUT      : po1    - The object. 
*              aepsge - Geometry resolution.
*
*
*
* INPUT/OUTPUT:
*              cmax   - The level value.
*
* OUTPUT     : pintdat - Maximum points  found. 
*              jstat  - status messages  
*                                = 1 : Maximumpoints found equal to level value
*                                = 2 : Maximumpoints found over level value
*                                = 0 : no maximum 
*                                < 0 : error
*
*
* METHOD     : This function is treating  point maximum. 
*              Otherwise it is computing edge/endpoint maximum  
*              by recurson on the edges. 
*              The maxima in the inner of the object is found by
*              calling another recursiv function.        
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              s1435      - Pick edge curve from a surface. 
*              s1438      - Pick endpoint from a curve.
*              s1162      - Find the intersections in the inner of the objects.
*              s1190      - A boxtest.
*              s6idnpt    - Put a new intpt to given intdat.
*              s6idput    - Put contence of one intdat in an other intdat.
*              s6idlis    - Compute list elements from given intdat.
*              s6idedg    - Uppdate an edge from given intdat.
*              freeEdge   - Free space occupied by given edge.
*              freePoint  - Free space occupied by given point.
*              freeCurve  - Free space occupied by given curve.
*              freeObject - Free space occupied by given object.
*              freeIntdat - Free space occupied by given intdat.
*              newEdge    - Create new edge-structure.
*              newIntpt   - Create new maximum point-structure.
*              newObject  - Create new object.
*
* WRITTEN BY : Ulf J. Krystad , 05.89.
*
*********************************************************************
*/                                     
{
  int    klevel=0;     /* Parameter into s1162                   */
  int    knum=0;       /* Parameter into s1162                   */
  int    kpar;         /* Fixed parameter direction.             */
  int    ki;           /* Counter.                               */    
  int    kedge;        /* Number of edges.                       */
  int idim  = 1;       /* Local dimension, always = 1            */
  int kstat = 0;       /* Local status variable.                 */
  int kpos  = 0;       /* Position of error.                     */
  double tpar;         /* Help variable used for parameter value
			  and geometric distance.                */
  SISLEdge   *qedge[2];        /* Edges for use in s1162().      */
  SISLObject *qdum = SISL_NULL;     /* Dummy  pointer.                */
  SISLObject *qob  = SISL_NULL;     /* Objects for use in recurson.   */
  SISLIntdat *qintdat = SISL_NULL;  /* Intdat for use in recurson.    */
  
  qedge[0] = SISL_NULL;
  qedge[1] = SISL_NULL;
  
  if (po1->iobj == SISLPOINT) 
    {
      /* It's a point, treat the case here and return. */
      
      /* Control the dimension. */
      if (po1->p1->idim != idim ) goto err106;
      
      /* Computing the distance beetween the point and level value. */
      tpar = po1->p1->ecoef[0] - *cmax;
      
      if (fabs(tpar) <= aepsge)
	
        /* The point is close enough to the level value to be a max. */
	*jstat = 1;         /* Mark maximum found. */
      
      else if (tpar > (double)0.0)
	{
	  
	  /* The point is greater than the level value . */
	  *jstat = 2;         /* Mark new maximum found. */
	  *cmax   = po1->p1->ecoef[0];
	}
      
      else 
	
	*jstat = 0;         /* Mark no maximum found. */
      
      
      if ( *jstat > 0 )
	{
	  SISLIntpt *qt;
	  
	  /* Add maximum  point. */
	  qt = newIntpt(0,cmax,DZERO);
	  if (qt == SISL_NULL) goto err101;
	  
	  /* Uppdate pintdat. */
	  s6idnpt(pintdat,&qt,1,&kstat);
	  if (kstat < 0) goto error;
	}
      
    }
  
  
  else if (po1->iobj > SISLPOINT)
    {
      /* It's a higher order geometry, treat the edges here and
	 use a recursiv function to treat the inner of the object       */
      
      
      *jstat = 0;
      /* Perform a boxtest */
      s1190(po1,cmax,aepsge,&kstat);
      if (kstat == 1) goto out;
      
      
      /*Create a dummy object, to be used when calling 
	the intersection routines
	treating two objects.*/
      if ((qdum = newObject(SISLPOINT)) == SISL_NULL) goto err101;
      
      
      
      kedge  = 2 * po1->iobj;
      kpar   = kedge/2;
      
      /* Create correct number of edges. */
      if ((qedge[0] = newEdge(kedge)) == SISL_NULL) goto err101;
      
      
      for (ki=0; ki<kedge; ki++)
	{
	  
	  /* Set  correct parameter direction to keep constant         */
	  kpar   = ((ki == kedge/2) ? kedge/2-1:kpar-1);
	  
	  /* Create one lower order helpobject */
	  if ((qob = newObject(po1->iobj - 1)) == SISL_NULL) goto err101;	
	  
	  
	  if (po1->iobj == SISLCURVE)
	    
	    /* Pick out end point from a curve. */
	    s1438(po1->c1,ki,&(qob->p1),&tpar,&kstat);
	  
	  else if (po1->iobj == SISLSURFACE)
	    
	    /* Pick out edge curve from a surface. */
	    s1435(po1->s1,ki,&(qob->c1),&tpar,&kstat);
	  
	  else
	    /* Unknown higher order object . */
	    goto err121;
	  
	  if (kstat < 0) goto error;
	  
	  /* Recursiv computing of end maximum. */
	  s1161(qob,cmax,aepsge,&qintdat,&kstat);
	  if (kstat < 0) goto error;
	  
	  if (kstat == 2)
	    {
	      
	      /* New maximum found, delete old ones */
	      if (*pintdat != SISL_NULL)
		{
		  freeIntdat(*pintdat);
		  *pintdat = SISL_NULL;
		}
	      
	      if (qedge[0] != SISL_NULL)
		{
		  /*  Empty the edges */
		  freeEdge(qedge[0]);
		  if ((qedge[0] = newEdge(kedge)) == SISL_NULL) goto err101;	      
		}
	      
	    }  
	  
	  
	  if (kstat)
	    {
	      /* Maximum found, add them to the list */
	      
	      *jstat = max(*jstat,kstat);         /* Mark maximum found. */
	      
	      /* Put maximum found on edges into pintdat. */
	      
	      /* Set parameter border values of object. */
	      s6idput(pintdat,qintdat,kpar,tpar,&kstat);
	      if (kstat < 0) goto error;
	      
	      /* Uppdate edge structure. */
	      s6idedg(po1,qdum,1,kpar+1,tpar,*pintdat,
		      &(qedge[0]->prpt[ki]),&(qedge[0]->ipoint),&kstat);
	      if (kstat < 0) goto error;
	    }
	  
	  if (qintdat != SISL_NULL) freeIntdat(qintdat);
	  qintdat = SISL_NULL;
	  freeObject(qob);
	}	  
      
      
      /* ---------------------------------------------------------------*/
      /* Treat the inner of higher order objects. */ 
      
      /* Before we enter internal maximum and subdivision we
	 initiate pointers to top level objects. */
      
      if (po1->o1 == SISL_NULL) po1->o1 = po1;
      
      /* Find the maximums in the inner of the object.  */
      s1162(po1,cmax,aepsge,pintdat,qedge,klevel,knum,&kstat);
      if (kstat < 0)  goto error;
      *jstat = max(*jstat,kstat);
      
      /* Organize the list in pintdat. */
      s6idlis(po1,po1,pintdat,&kstat);
      if (kstat < 0)  goto error;
    }
  
  else 
    /* Unknown  object . */
    goto err121;
  
  
  goto out; 
  
  
  
  /* -------------- ERROR HANDLING ----------------------------------------*/
  
  /* Error in space allocation.  */
 err101: *jstat = -101;
  s6err("s1161",*jstat,kpos);
  goto out;
  
  /* Error. Dimensions conflicting.  */
 err106: *jstat = -106;
  s6err("s1161",*jstat,kpos);
  goto out;
  
  /* Error. Kind of object does not exist.  */
 err121: *jstat = -121;
  s6err("s1161",*jstat,kpos);
  goto out;
  
  /* Error in lower order routine.  */
  error : *jstat = kstat;
  s6err("s1161",*jstat,kpos);
  goto out;
  
 out:
  /* Free the edges used in s1162. */
  if (qedge[0] != SISL_NULL) freeEdge(qedge[0]);
  
  /* Free the dummy object(point). */
  if (qdum != SISL_NULL) freeObject(qdum);
  
}

