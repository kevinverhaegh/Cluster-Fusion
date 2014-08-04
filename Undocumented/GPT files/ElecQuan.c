/* ElecQuan.c - Creates an electron particle distribution, based on the ion particle distribution */

#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <iostream>
#include <random>
#include "elem.h"
#include <time.h> 

extern double dblpulsarrand(void) ;

void ElecQuan_init(gptinit *init)
{
  gptparset *set1 ;
  gptinitpar *par1 ;
  gptparset *set2 ;
  gptinitpar *par2 ;
  double Rmin, Rmax ;
  double Eeig, Rxd, Ryd, Rzd, Rx, Ry, Rz, Epot, Ekin, vR, r ;
  char *name1 ;
  char *name2 ;
  int i, len ;

  Eeig = 13.6 * gpt_qe;

  std::tr1::mt19937 eng;
  eng.seed((unsigned int)time(NULL));
  
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(set1,set2, Rmin, Rmax)\n", gptgetname(init) ) ;

  name1 = gptgetargstring(init,1) ;
  name2 = gptgetargstring(init,2) ;
  Rmin = gptgetargdouble(init,3) ;
  Rmax = gptgetargdouble(init,4) ;

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
	r = (2*dblpprand()-1)*(Rmax - Rmin) + Rmin;
	printf("random: %E \n", 2*dblpprand()-1) ; 
	printf("r: %E \n", r);
	printf("Rmax: %E \n", Rmax);
	printf("Rmin: %E \n", Rmin);
	
    	std::normal_distribution<double> distribution(0.0,1.0);
	
	Rxd = distribution(eng);
	Ryd = distribution(eng);
	Rzd = distribution(eng);
	printf("Rxd, Ryd, Rzd: %E, %E, %E \n", Rxd, Ryd, Rzd);
	
	Rx = (Rxd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd))*r ;
	Ry = (Ryd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd))*r ;
	Rz = (Rzd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd))*r ;
	printf("Rx, Ry, Rz: %E, %E, %E \n", Rx, Ry, Rz);
	
    par2[i].Wr[0] = par1[i].Wr[0] + Rx;
    par2[i].Wr[1] = par1[i].Wr[1] + Ry;
    par2[i].Wr[2] = par1[i].Wr[2] + Rz;
	
	Epot = (gpt_qe*gpt_qe / (4*gpt_pi*gpt_eps0))*(1/r) ;
	printf("Epot: %E \n", Epot) ;
	
	Ekin = Eeig + Epot;
	printf("Ekin: %E \n", Ekin) ;
	vR = sqrt(2*Ekin / gpt_me) ;
	printf("vR: %E \n", vR) ;
	par2[i].GBr[0] = (vR/gpt_c)*(Rxd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd)) ;
	par2[i].GBr[1] = (vR/gpt_c)*(Ryd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd)) ;
	par2[i].GBr[2] = (vR/gpt_c)*(Rzd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd)) ;

	printf("vx: %E \n", vR*(Rxd/sqrt(Rxd*Rxd + Ryd*Ryd + Rzd*Rzd))) ;
  }
}
