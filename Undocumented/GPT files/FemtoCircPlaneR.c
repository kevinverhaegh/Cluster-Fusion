/* FemtoCircPlaneR.c: Femtosecond pulse circular polarization righthand */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct FemtoCircPlaneR_info
{
  double E0 ;
  double omega ;
  double tau ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int FemtoCircPlaneR_sim(gptpar *par,double t,struct FemtoCircPlaneR_info *info) ;

/* Initialization routine */
void FemtoCircPlaneR_init(gptinit *init)
{
  struct FemtoCircPlaneR_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=3 )
    gpterror( "Syntax: %s(ECS,E0,omega,tau)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct FemtoCircPlaneR_info *)gptmalloc( sizeof(struct FemtoCircPlaneR_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->E0       = gptgetargdouble(init,1) ;
  info->omega    = gptgetargdouble(init,2) ;
  info->tau      = gptgetargdouble(init,3) ;

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  gptaddEBelement( init, FemtoCircPlaneR_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int FemtoCircPlaneR_sim(gptpar *par,double t,struct FemtoCircPlaneR_info *info)
{
  /* Copy of parameters in info structure for convenience */
  double E0, omega, tau ;

  /* Retrieve parameters from info structure */
  E0       = info->E0 ;
  omega    = info->omega ;
  tau      = info->tau ;

  /* Calculate electromagnetic fields from the above parameters
   * Particle coordinates: X,Y and Z   must be written UPPERCASE
   * Simulation time     : t           must be written lowercase
   * Electric field      : EX = ... ; EY = ...  ; EZ = ... ;
   * Magnetic field      : BX = ... ; BY = .... ; BZ = ... ;
   */
  EX = (1/sqrt(2.0))*E0*exp(-(1/(2*tau*tau))*(t - Z/gpt_c))*sin(omega*(t - Z/gpt_c)) ;
  EY = (1/sqrt(2.0))*E0*exp(-(1/(2*tau*tau))*(t - Z/gpt_c))*cos(omega*(t - Z/gpt_c)) ;
  BX = -(1/(gpt_c*(sqrt(2.0))))*E0*exp(-(1/(2*tau*tau))*(t - Z/gpt_c))*cos(omega*(t - Z/gpt_c)) ;
  BY = (1/(gpt_c*(sqrt(2.0))))*E0*exp(-(1/(2*tau*tau))*(t - Z/gpt_c))*sin(omega*(t - Z/gpt_c)) ;

  /* Return 1 to notify particle is INSIDE element */
  return( 1 ) ;
}
