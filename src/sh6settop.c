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
 * $Id: sh6settop.c,v 1.2 2001-03-19 16:06:03 afr Exp $
 *
 */


#define SH6SETTOP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
      sh6settop(SISLIntpt *pt,int ilist,int left1,int right1,int left2,int right2,int *jstat)
#else
void sh6settop(pt,ilist,left1,right1,left2,right2,jstat)
   SISLIntpt *pt;
   int       ilist;
   int       left1;
   int       right1;
   int       left2;
   int       right2;
   int       *jstat;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Given an Intpt and a list which it lies in, set the
*              pre-topology information. If list does not
*              exist give an error message.
*
*
* INPUT      : pt       - Pointer to the Intpt.
*              ilist    - Index specifying a list.
*              left1    - pre-topology data.
*              right1   - pre-topology data.
*              left2    - pre-topology data.
*              right2   - pre-topology data.
*
*
* OUTPUT     : pt       - Pointer to the updated Intpt.
*              jstat    - Error flag.
*                         jstat =  0  => OK.
*                         jstat = -1  => ilist is out of range.
*                         jstat = -2  => Error.
*
*
* METHOD     : 
*
*
* REFERENCES :
*
* WRITTEN BY : Michael Floater, SI, Oslo, Norway. June 91.
*
*********************************************************************
*/
{
   *jstat=0;

   /* Check pt. */

   if(pt == SISL_NULL) goto err2;

   /* Check ilist. */

   if(ilist >= 0 && ilist < pt->no_of_curves)
   {
       pt->left_obj_1[ilist]=left1;
       pt->right_obj_1[ilist]=right1;
       pt->left_obj_2[ilist]=left2;
       pt->right_obj_2[ilist]=right2;
   }
   else if(pt->no_of_curves == 0 && ilist == 0)
   {
       pt->left_obj_1[0]=left1;
       pt->right_obj_1[0]=right1;
       pt->left_obj_2[0]=left2;
       pt->right_obj_2[0]=right2;
   }
   else if(ilist == -1)
   {
       pt->left_obj_1[0]=left1;
       pt->right_obj_1[0]=right1;
       pt->left_obj_2[0]=left2;
       pt->right_obj_2[0]=right2;
   }
   else goto err1;


   /* Data is set. */

   goto out;
   

err1:
   /* Error. ilist is out of range. */
   
   *jstat = -1;
   s6err("sh6settop",*jstat,0);
   goto out;

err2:
   /* Error in input. pt is SISL_NULL. */
   
   *jstat = -2;
   s6err("sh6settop",*jstat,0);
   goto out;
   
   
   out :
      return;
}
