/* translateset.c - Copies particle set and translates it */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

extern double dblpulsarrand(void) ;

void translateset_init(gptinit *init)
{
  double x0,y0,z0 ;
  gptparset *set1 ;
  gptinitpar *par1 ;
  gptparset *set2 ;
  gptinitpar *par2 ;
  char *name1 ;
  char *name2 ;
  int i, len ;

  if( gptgetargnum(init)!=5 )
    gpterror( "Syntax: %s(set1,set2,x0,y0,z0)\n", gptgetname(init) ) ;

  name1 = gptgetargstring(init,1) ;
  name2 = gptgetargstring(init,2) ;
  x0    = gptgetargdouble(init,3) ;
  y0    = gptgetargdouble(init,4) ;
  z0    = gptgetargdouble(init,5) ;

  /* Get particle set 1 */
  if( gpttestparset( name1 )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name1 ) ;
  set1 = gptgetparset( name1 ) ;
  par1 = gptgetparsetpars( set1,&len ) ;
  
   /* Get particle set 2 */
  if( gpttestparset( name2 )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name2 ) ;
  set2 = gptgetparset( name2 ) ;
  par2 = gptgetparsetpars( set2,&len ) ;

  /* Copies set 1 and puts it in set 2 with translation */
  for( i=0 ; i<len ; i++ )
  {

    par2[i].Wr[0] = par1[i].Wr[0] +x0;
    par2[i].Wr[1] = par1[i].Wr[1] +y0 ;
    par2[i].Wr[2] = par1[i].Wr[2] +z0;
  }
}
