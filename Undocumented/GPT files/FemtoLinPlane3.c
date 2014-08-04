/* FemtoLinPlane3.c:  */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct FemtoLinPlane3_info
{
  double wd ;
  double tau ;
  double E0 ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int FemtoLinPlane3_sim(gptpar *par,double t,struct FemtoLinPlane3_info *info) ;

/* Initialization routine */
void FemtoLinPlane3_init(gptinit *init)
{
  struct FemtoLinPlane3_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,wd,tau,E0)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct FemtoLinPlane3_info *)gptmalloc( sizeof(struct FemtoLinPlane3_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->wd       = gptgetargdouble(init,1) ;
  info->tau      = gptgetargdouble(init,2) ;
  info->E0       = gptgetargdouble(init,3) ;

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  gptaddEBelement( init, FemtoLinPlane3_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int FemtoLinPlane3_sim(gptpar *par,double t,struct FemtoLinPlane3_info *info)
{
  /* Copy of parameters in info structure for convenience */
  double wd, tau, E0;

  /* Retrieve parameters from info structure */
  wd       = info->wd ;
  tau      = info->tau ;
  E0       = info->E0 ;

  /* Calculate electromagnetic fields from the above parameters
   * Particle coordinates: X,Y and Z   must be written UPPERCASE
   * Simulation time     : t           must be written lowercase
   * Electric field      : EX = ... ; EY = ...  ; EZ = ... ;
   * Magnetic field      : BX = ... ; BY = .... ; BZ = ... ;
   */

  EX = E0 * exp(-(1/(2*(tau*tau)))*(t- Z/gpt_c- 9.8735e-14)*(t-Z/gpt_c - 9.8735e-14))*cos(wd*(t- Z/gpt_c-9.8735e-14)) ;
  BY = (E0/gpt_c) * exp(-(1/(2*(tau*tau)))*(t- Z/gpt_c- 9.8735e-14)*(t-Z/gpt_c- 9.8735e-14))*cos(wd*(t- Z/gpt_c-9.8735e-14)) ;

  /* Return 1 to notify particle is INSIDE element */
  return( 1 ) ;
}
