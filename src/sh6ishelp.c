/*****************************************************************************/
/*                                                                           */
/*                                                                           */
/* (c) Copyright 1989,1990,1991,1992 by                                      */
/*     Senter for Industriforskning, Oslo, Norway                            */
/*     All rights reserved. See the copyright.h for more details.            */
/*                                                                           */
/*****************************************************************************/

#include "copyright.h"

/*
 *
 * $Id: sh6ishelp.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define SH6ISHELP

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
int 
      sh6ishelp(SISLIntpt *pt)
#else
int sh6ishelp(pt)
    SISLIntpt *pt;
#endif   
/*
*********************************************************************
*
*********************************************************************
*
* PURPOSE    : Boolean. Is pt a help point?
*
*
* INPUT      : pt       - Pointer to the Intpt.
*
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
   int flag = 0;

   if(pt != NULL && pt->iinter < 0) flag = 1;


   goto out;
   

   
   out :
      return flag;
}






