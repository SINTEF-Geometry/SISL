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
 * $Id: s6line.c,v 1.1 1994-04-21 12:10:42 boh Exp $
 *
 */


#define S6LINE

#include "sislP.h"

#if defined(SISLNEEDPROTOTYPES)
void 
s6line(double epoint[])
#else
void s6line(epoint)
     double epoint[];
#endif
{
  printf("\n s6line: %f    %f    %f ",epoint[0],epoint[1],epoint[2]);
}
