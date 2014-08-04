/* setellipse2.c - Set homogeneous ellipse */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include "elem.h"

extern double dblpulsarrand(void) ;

void setellipse2_init(gptinit *init)
{
  double a,b,c ;
  double x,y,z ;
  double x0,y0,z0 ;
  gptparset *set ;
  gptinitpar *par ;
  char *name ;
  int i, len ;

  if( gptgetargnum(init)!=7 )
    gpterror( "Syntax: %s(set,a,b,c)\n", gptgetname(init) ) ;

  name = gptgetargstring(init,1) ;
  a    = gptgetargdouble(init,2) ;
  b    = gptgetargdouble(init,3) ;
  c    = gptgetargdouble(init,4) ;
  x0   = gptgetargdouble(init,5) ;
  y0   = gptgetargdouble(init,6) ;
  z0   = gptgetargdouble(init,7) ;

  /* Get particle set */
  if( gpttestparset( name )==NULL )
    gptwarning( "The particle set \"%s\" does not exist\n", name ) ;
  set = gptgetparset( name ) ;
  par = gptgetparsetpars( set,&len ) ;

  /* Set ellipse */
  for( i=0 ; i<len ; i++ )
  {
    do
    {
      /* Uniform in box between -a and a, -b and b, -c and c */
      x = 2*dblpprand()-1 ;
      y = 2*dblpprand()-1 ;
      z = 2*dblpprand()-1 ;
    } 
    while( x*x+y*y+z*z >= 1 ) ;

    par[i].Wr[0] = a*x +x0 ;
    par[i].Wr[1] = b*y +y0 ;
    par[i].Wr[2] = c*z +z0 ;
  }
}
