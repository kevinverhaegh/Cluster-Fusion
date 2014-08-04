/* ADKIonize.c:  */
/* ADKIonize.c:  */

#include <stdio.h>
#include <math.h>
#include "elem.h"
#include <cstdlib>

/* Info structure containing all relevant parameters for this element */
struct ADKIonize_info 
{
  gptparset *set ;
  struct axis *axis ;
  double rmacro ;
  int Ato ;
  int Nr ;
  double Ez[200];
} ;

/* Forward declaration of the routine calculating the electromagnetic fields */
static int ADKIonize_out(double t, double *dt, double *x, void *vinfo ) ;

/* Initialization routine */
void ADKIonize_init(gptinit *init)
{
  struct ADKIonize_info *info ;

  /* Read Element Coordinate System (ECS) from parameter list */
  gptbuildECS( init ) ;

  /* Print usage line when the number of parameters is incorrect */
  if( gptgetargnum(init)!=4 )
    gpterror( "Syntax: %s(ECS,set,rmacro,Ato,Nr)\n", gptgetname(init) ) ;

  /* Allocate memory for info structure */
  info = (struct ADKIonize_info *)gptmalloc( sizeof(struct ADKIonize_info) ) ;

  /* Read all parameters as doubles and store them in info structure */
  info->axis   = getaxis("wcs") ;
  info->set    = gptgetparset(gptgetargstring(init,1)) ;
  info->rmacro = gptgetargdouble(init,2) ;
  info->Ato    = gptgetargint(init,3) ;
  info->Nr     = gptgetargint(init,4) ;

 FILE *infile;
 char str[9999], Ezch[7] ;

printf("Reading ionization data \n");

  for(int i=0; i<info->Ato; i++)
  {
	// Read ionization energies
	infile = fopen("Ez.txt", "r");

	if (infile == NULL) perror ("Error opening file") ; 
	
	//Get string of atom needed in ionization energy file

	if (info->Ato != 1)
	{
		for (int j=1 ; j<info->Nr ; j++)
		{
		fgets(str, 9999, infile);
		}
	}
   

	
	//Get following ionization energy (e.q. get ionization energy of the Z+1 level)

	fgets(str, 8*i, infile) ;
	fgets(Ezch, 8, infile) ;
	double Ezv = atof(Ezch);
	fclose(infile) ;
	
	info->Ez[i]=Ezv ; 
	printf("Ionization stage found \n");
	 
 }

  /* Register the routine calculating the electromagnetic fields to the GPT kernel */
  odeaddoutfunction( ODEFNC_USR, ADKIonize_out, info ) ;
}


/* The following routine calculates the electromagnetic fields */
static int ADKIonize_out(double t, double *dt, double *x, void *vinfo )
{
  struct ADKIonize_info *info = (struct ADKIonize_info *)vinfo ;
  
  // Atomic physics parameters
  double auE = 5.14220652e11;
  double auEn = 4.35974417e-18;
  double auL = 5.2917720859e-11;
  double aut = 2.418884326505e-17;

  //printf("Starting ADKIonize \n");

  // Loop over all particles
  for(int i=0 ; i<numpar ; i++)
    if( (pars[i].alive) && (pars[i].q>=0 ))
  {

	// Determine Z 
	double qi = pars[i].q ;
	int Zi = int ( - qi / gpt_qe);
	double Ez=info->Ez[Zi] ; 
	//double Ez=14;
	if (info->Ato != Zi) 
	{
	//printf("Z, %f \n",Zi);

	//Printf debugging\
	
	//printf("Ionizable element found \n");
	//printf("Field x component:, %e \n", pars[i].WE[0]);
	//printf("x position:, %e \n", pars[i].Wr[0]);
	//printf("charge:, %e \n", pars[i].q);
	//printf("mass:, %e \n", pars[i].m);
	//printf("ID:, %e \n", pars[i].ID);
    	//printf("Time:,%e \n", t);
	//printf("dt: %e \n", *dt);
		
	// Calculate electric field strength
	
	double Ef = sqrt(pars[i].WE[0]*pars[i].WE[0] + pars[i].WE[1]*pars[i].WE[1] + pars[i].WE[2]*pars[i].WE[2]) ;
	//printf("Field strength, %e \n", Ef); 

	//printf("Ionization energy %f \n", Ez);	
	//double Ez=14;
	//Ionization energy of next level now known

	//Check if ionization occurs with ADK model

	//Calculate ADK tunnel probability

	//Conversion to atomic units

	double EzA = (-gpt_qe * Ez)/auEn;
	double EfA = Ef/auE; 

	//Calculate paramters in ADK formula
	
	double nstar = (1.0+ Zi)/sqrt(2*EzA);
	//printf("nstar: %e\n", nstar);

	double Cnstar=pow(((2.0* exp(1.0))/nstar), nstar)/sqrt(2*gpt_pi*nstar);
	//printf("Cnstar: %e\n", Cnstar);
	
	//Calculating factorials
	
	int c;
	
	//assume l=m=0
	int l=0;
	int m=0;
	
	double fac1=1;
	double fac2=1;
	double fac3=1;
	
	for (c = 1; c <= l+abs(m); c++)
	{
		fac1=fac1*c;
	}
	
	for (c = 1; c <= l - abs(m); c++)
	{
		fac2=fac2*c;
	}
	
	for (c = 1; c <= abs(m); c++)
	{
		fac3=fac3*c;
	}
	
	double flm =  ((2.0*l+1)*fac1)/(pow(2.0, abs(m)))*fac3*fac2;
	//printf("flm: %e\n", flm);
	
	double w = Cnstar*Cnstar*flm*(((1.0+Zi)*(Zi+1.0))/(2.0*nstar*nstar))*sqrt((nstar*nstar*nstar)/((Zi+1.0)*(Zi+1.0)*(Zi+1.0)))*pow((2.0*(Zi+1.0)*(Zi+1.0)*(Zi+1.0))/(EfA*nstar*nstar*nstar),(2*nstar-abs(m)-1))*exp(-(2*(Zi+1.0)*(Zi+1.0)*(Zi+1.0))/(3*nstar*nstar*nstar*EfA));
	//printf("dw/dt: %e\n", w);	

	double tstep = *dt;
	//printf("dt: %e\n", tstep);
	
	double proba=w*(tstep/aut);
	//printf("Ionization probability:,%e \n", proba);
	
	double randomV = dblpprand();
	//printf("Random number:,%e \n", randomV);
	
		if (proba>randomV)
		{
	
		//Ionization occurs
		//printf("Ionization occurs \n"); 
	
		//Add charge to the ion
		pars[i].q = pars[i].q - gpt_qe ;
	
		//create electron
		double Wr[3], GBr[3] ;

		//Electron position according to ADK top

		double displ = - auL * sqrt((Zi+1)/(Ef/auE)); 
	
		Wr[0] = pars[i].Wr[0] + displ * (pars[i].WE[0] / sqrt(pars[i].WE[0] * pars[i].WE[0] + pars[i].WE[1] * pars[i].WE[1] + pars[i].WE[2] *pars[i].WE[2]));
		Wr[1] = pars[i].Wr[1] + displ * (pars[i].WE[1] / sqrt(pars[i].WE[0] * pars[i].WE[0] + pars[i].WE[1] * pars[i].WE[1] + pars[i].WE[2] *pars[i].WE[2]));
		Wr[2] = pars[i].Wr[2] + displ * (pars[i].WE[2] / sqrt(pars[i].WE[0] * pars[i].WE[0] + pars[i].WE[1] * pars[i].WE[1] + pars[i].WE[2] *pars[i].WE[2]));
		
		//Electron in rest
		
		GBr[0] = 0;
		GBr[1] = 0;
		GBr[2] = 0;

		gptsimaddparmqnartid(info->set, Wr, GBr, gpt_me, gpt_qe, 1, info->axis,info->rmacro, t, 0) ;
		}

	}	

	
	

 }
  return 0 ;
}
