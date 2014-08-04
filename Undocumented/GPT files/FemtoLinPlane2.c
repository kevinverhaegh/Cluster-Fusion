/* FemtoLinPlane2.c:  */

#include <stdio.h>
#include <math.h>
#include "elem.h"

/* Info structure containing all relevant parameters for this element */
struct FemtoLinPlane2_info
{
  double wd ;
  double tau ;
  double E0 ;
  double ZL ;
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int FemtoLinPlane2_sim(gptpar *par,double t,struct FemtoLinPlane2_info *info) ;

/* Initialization routine */
void FemtoLinPlane2_init(gptinit *init)
{
  struct FemtoLinPlane2_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,wd,tau,E0,ZL)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct FemtoLinPlane2_info *)gptmalloc( sizeof(struct FemtoLinPlane2_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->wd       = gptgetargdouble(init,1) ;
  info->tau      = gptgetargdouble(init,2) ;
  info->E0       = gptgetargdouble(init,3) ;
  info->ZL       = gptgetargdouble(init,4) ;

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  gptaddEBelement( init, FemtoLinPlane2_sim, gptfree, GPTELEM_LOCAL, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int FemtoLinPlane2_sim(gptpar *par,double t,struct FemtoLinPlane2_info *info)
{
  /* Copy of parameters in info structure for convenience */
  double wd, tau, E0, ZL ;

  /* Retrieve parameters from info structure */
  wd       = info->wd ;
  tau      = info->tau ;
  E0       = info->E0 ;
  ZL       = info->ZL ;

  /* Calculate electromagnetic fields from the above parameters
   * Particle coordinates: X,Y and Z   must be written UPPERCASE
   * Simulation time     : t           must be written lowercase
   * Electric field      : EX = ... ; EY = ...  ; EZ = ... ;
   * Magnetic field      : BX = ... ; BY = .... ; BZ = ... ;
   */

  EY = E0 * exp(-(1/(2*(tau*tau)))*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c))*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c)))*cos(wd*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c))) ;
  BX = (E0/gpt_c) * exp(-(1/(2*(tau*tau)))*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c))*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c)))*cos(wd*(t-2*(2.0*sqrt(log(2.0)))*tau - (ZL/gpt_c + Z/gpt_c))) ;

  /* Return 1 to notify particle is INSIDE element */
  return( 1 ) ;
}
