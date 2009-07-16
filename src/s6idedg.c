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
 * $Id: s6idedg.c,v 1.2 2001-03-19 15:59:01 afr Exp $
 *
 */


#define S6IDEDG

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6idedg(SISLObject *po1,SISLObject *po2,int iobj,int ipar,double apar,
	     SISLIntdat *pintdat,SISLPtedge **rptedge,int *jnum,int *jstat)
#else
void s6idedg(po1,po2,iobj,ipar,apar,pintdat,rptedge,jnum,jstat)
     SISLObject *po1;
     SISLObject *po2;
     int    iobj;
     int    ipar;
     double apar;
     SISLIntdat *pintdat;
     SISLPtedge **rptedge;
     int    *jnum;
     int    *jstat;
#endif
/*
*********************************************************************
*
*********************************************************************
*                                                                   
* PURPOSE    : To make a list of ptedges pointing to intersection
*              points with one constant parameter value apar and if it
*              exist more than one parameter the other value have to
*              be between tstart and tend. A pointer to this list
*              is returned by rptedge.
*
*
*
* INPUT      : po1      - First object in intersection.
*              po2      - Second object in intersection.
*              iobj     - Number of object to pick edge of.
*              ipar     - Number of parameter direction to pick edge  of.
*              apar     - The edge/end parameter value.
*              pintdat  - Pointer to intersection data.
*
*
* OUTPUT     : rptedge  - Pointer to a pointer to the ptedge list.
*              jnum     - += number of elements in list.
*              jstat    - status messages  
*                               = 0      : OK!
*                               < 0      : error
*
*
* METHOD     : 
*
*
* REFERENCES :
*
*-
* CALLS      : s6err      - Gives error message.
*              s6idnpt    - Insert a new intpt structure.
*              copyIntpt  - Copy an intpt structure.
*              newIntdat  - Create new intdat structure.
*
* WRITTEN BY : Arne Laksaa, 05.89.
* Revised by : Christophe Rene Birkeland, SINTEF Oslo, May 1993.
*              Error message corrected
*
*********************************************************************
*/                                     
{
  int kpos=0;                /* Position of error.                       */
  int kpar=0;                /* Numper of parameter direction second obj.*/
  int ki,kj;                 /* Counters                                 */
  double sstart[4],send[4];  /* Parameter boarders on the other obj.     */
  SISLPtedge *pte = SISL_NULL;    /* Pointers to new ptedge.                  */
  
  /* Initiate to emty list. */
  
  *rptedge = SISL_NULL;  
  *jstat = 0;
  
  /* We have to be sure that we have an intdat structure. */
  
  if (pintdat == SISL_NULL) goto out;
  
  /* Uppdate parameter boarder. */
  
  if (po1->iobj == SISLCURVE)
    {
      if (iobj == 1)
        {
	  sstart[0] = apar;
	  send[0]   = apar;
        }
      else
        {
	  sstart[0] = po1->c1->et[po1->c1->ik - 1];
	  send[0]   = po1->c1->et[po1->c1->in];
        }
      kpar = 1;
    }
  else if (po1->iobj == SISLSURFACE)
    {
      if (iobj == 1 && ipar == 1)
        {
	  sstart[0] = apar;
	  send[0]   = apar;
        }
      else
        {
	  sstart[0] = po1->s1->et1[po1->s1->ik1 - 1];
	  send[0]   = po1->s1->et1[po1->s1->in1];
        }
      if (iobj == 1 && ipar == 2)
        {
	  sstart[1] = apar;
	  send[1]   = apar;
        }
      else
        {
	  sstart[1] = po1->s1->et2[po1->s1->ik2 - 1];
	  send[1] = po1->s1->et2[po1->s1->in2];
        }
      kpar = 2;
    }
  
  
  if (po2->iobj == SISLCURVE)
    {
      if (iobj == 2)
        {
	  sstart[kpar] = apar;
	  send[kpar]   = apar;
        }
      else
        {
	  sstart[kpar] = po2->c1->et[po2->c1->ik - 1];
	  send[kpar]   = po2->c1->et[po2->c1->in];
        }
    }
  else if (po2->iobj == SISLSURFACE)
    {
      if (iobj == 2 && ipar == 1)
        {
	  sstart[kpar] = apar;
	  send[kpar]   = apar;
        }
      else
        {
	  sstart[kpar] = po2->s1->et1[po2->s1->ik1 - 1];
	  send[kpar] = po2->s1->et1[po2->s1->in1];
        }
      if (iobj == 2 && ipar == 2)
        {
	  sstart[kpar+1] = apar;
	  send[kpar+1]   = apar;
        }
      else
        {
	  sstart[kpar+1] = po2->s1->et2[po2->s1->ik2 - 1];
	  send[kpar+1]   = po2->s1->et2[po2->s1->in2];
        }
    }
    
  /* We have to go trough all intersection points to search for edges. */
  
  for (ki=0; ki<pintdat->ipoint; ki++)
    {
      for (kj=0; kj<pintdat->vpoint[ki]->ipar; kj++)
        if ((DEQUAL(sstart[kj],pintdat->vpoint[ki]->epar[kj]) ||
	     sstart[kj] < pintdat->vpoint[ki]->epar[kj]) &&
	    (DEQUAL(send[kj],pintdat->vpoint[ki]->epar[kj]) ||
	     send[kj] > pintdat->vpoint[ki]->epar[kj]));
	else
	  goto end;
      
      if (pte == SISL_NULL)
        {
	  pte = newPtedge(pintdat->vpoint[ki]);
	  if (pte == SISL_NULL) goto err101;
	  
	  (*rptedge) = pte;
	  
	  (*jnum)++;
        }
      else
        {
	  pte->pnext = newPtedge(pintdat->vpoint[ki]);
	  if (pte->pnext == SISL_NULL) goto err101;
	  
	  pte = pte->pnext;
	  
	  (*jnum)++;
        }
    end:;
    }
  
  goto out;
  
  /* Error in space allocation.  */

  err101: 
    *jstat = -101;
    s6err("s6idedg",*jstat,kpos);
    goto out;

 out:  ;
}
