/* EOscillatingSin.c:  */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct EOscillatingSin_info
{
  double E0 ;
  double wd ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int EOscillatingSin_sim(gptpar *par,double t,struct EOscillatingSin_info *info) ;

/* Initialization routine */
void EOscillatingSin_init(gptinit *init)
{
  struct EOscillatingSin_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=2 )
    gpterror( "Syntax: %s(ECS,E0,wd)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct EOscillatingSin_info *)gptmalloc( sizeof(struct EOscillatingSin_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->E0       = gptgetargdouble(init,1) ;
  info->wd       = gptgetargdouble(init,2) ;

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  gptaddEBelement( init, EOscillatingSin_sim, gptfree, GPTELEM_GLOBAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int EOscillatingSin_sim(gptpar *par,double t,struct EOscillatingSin_info *info)
{
  /* Copy of parameters in info structure for convenience */
  double E0, wd ;

  /* Retrieve parameters from info structure */
  E0       = info->E0 ;
  wd       = info->wd ;

  /* Calculate electromagnetic fields from the above parameters
   * Particle coordinates: X,Y and Z   must be written UPPERCASE
   * Simulation time     : t           must be written lowercase
   * Electric field      : EX = ... ; EY = ...  ; EZ = ... ;
   * Magnetic field      : BX = ... ; BY = .... ; BZ = ... ;
   */
  EY = E0*sin(wd*t) ;
  BX = (E0/gpt_c)*sin(wd*t) ;

  /* Return 1 to notify particle is INSIDE element */
  return( 1 ) ;
}
