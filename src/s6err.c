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
 * $Id: s6err.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6ERR

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void s6err(char *rut,int jstat,int ipos)
#else
void s6err(rut,jstat,ipos)
     char *rut;
     int  jstat;
     int  ipos;
#endif
{
   (void)fprintf(stderr,"\nError status : %d",jstat);
   (void)fprintf(stderr,"   Call from routine : %s",rut);
   (void)fprintf(stderr,"   Position : %d\n",ipos);
}
