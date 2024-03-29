#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "allvars.h"
#include "proto.h"

void gas_inflow(int p, double time)
{
  double r_in, r_out, gas_old[RNUM], newarea, r1, r2, frac, alpha, inflowfrac, rgas;
  int j, index;
	
#ifdef DETAILED_METALS_AND_MASS_RETURN
  struct metals gasmetal_old[RNUM];
#ifdef INDIVIDUAL_ELEMENTS
  int kk;
  double gasElements_old[RNUM][NUM_ELEMENTS];
#endif
#else
  double gasmetal_old[RNUM];
#endif

  if(Gal[p].ColdGas<1.0e-6)
    rgas=Gal[p].ColdGasRadius;
  else
    {
      rgas=0.5*RingRadius[0]*Gal[p].ColdGasRings[0];
      for(j=1;j<RNUM;j++) rgas+=(0.5*(RingRadius[j-1]+RingRadius[j])*Gal[p].ColdGasRings[j]);
      rgas=rgas/Gal[p].ColdGas/2.0;      //2.0=mean radius/scale length for exponential disk
    }
   //rgas=Gal[p].ColdGasRadius/3.0;
	 
  //inflowfrac=1.00;
  //alpha=1-1.34e3*time;
  for (j=0; j<RNUM; j++)
    {
      inflowfrac=1.00;
      gas_old[j]=Gal[p].ColdGasRings[j]*inflowfrac;
      gasmetal_old[j]=metals_add(metals_init(),Gal[p].MetalsColdGasRings[j],inflowfrac);
#ifdef INDIVIDUAL_ELEMENTS
      for (kk=0;kk<NUM_ELEMENTS;kk++)
	gasElements_old[j][kk]=Gal[p].ColdGasRings_elements[j][kk]*inflowfrac;
#endif

      Gal[p].ColdGasRings[j]*=(1-inflowfrac);
      Gal[p].MetalsColdGasRings[j]=metals_add(metals_init(),Gal[p].MetalsColdGasRings[j],(1.-inflowfrac));
#ifdef INDIVIDUAL_ELEMENTS
      for (kk=0;kk<NUM_ELEMENTS;kk++)
	Gal[p].ColdGasRings_elements[j][kk]=Gal[p].ColdGasRings_elements[j][kk]*(1.-inflowfrac);
#endif
    }

  for (j=0; j<RNUM; j++)
    {
      alpha=1-GasInflowVel*time/Hubble_h; //time unit: (Mpc/h)/(km/s)
      r_out=RingRadius[j]*alpha;
      if(j==0) r_in=0.0;
      else r_in=RingRadius[j-1]*alpha;
      if(r_in<1.0e-8) r_in=1.0e-8;

      //constant inflow velocity
      //vgas=3.0; // km/s
      //r_out=RingRadius[j]-vgas*time; //time unit: (Mpc/h)/(km/s); radius unit: Mpc/h
      //if(j==0) r_in=0.0;
      //else r_in=RingRadius[j-1]-vgas*time;

      //constant inflow velocity
      if(r_out<=RingRadius[0])
	{
	  index=0;
	  Gal[p].ColdGasRings[0]+=gas_old[j];
	  Gal[p].MetalsColdGasRings[0]=metals_add(Gal[p].MetalsColdGasRings[0],gasmetal_old[j],1.);
#ifdef INDIVIDUAL_ELEMENTS
	  for (kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].ColdGasRings_elements[0][kk]+=gasElements_old[j][kk];
#endif
	}
      else
	{
	  newarea=r_out*r_out-r_in*r_in;
	  //index=floor(log(1000*r_out)/log(1.2)+3); //The inverse of RingRadius[i]= pow(1.2,i-3)/1000
	  for(index=0;RingRadius[index]<=r_out&&index<RNUM-1;index++);	//ring number of the r_out
	  index--;
	  r1=RingRadius[index]; r2=r_out;
	  while(r1>=r_in)
	    {
	      frac=(r2*r2-r1*r1)/newarea;
	      Gal[p].ColdGasRings[index+1]+=frac*gas_old[j];
	      Gal[p].MetalsColdGasRings[index+1]=metals_add(Gal[p].MetalsColdGasRings[index+1],gasmetal_old[j],frac);
#ifdef INDIVIDUAL_ELEMENTS
	      for (kk=0;kk<NUM_ELEMENTS;kk++)
		Gal[p].ColdGasRings_elements[index+1][kk]+=(gasElements_old[j][kk]*frac);
#endif
	      index--;
	      r2=r1;
	      if(index>-1) r1=RingRadius[index];
	      else r1=0.0;
	    }
	  r1=r_in;
	  frac=(r2*r2-r1*r1)/newarea;
	  Gal[p].ColdGasRings[index+1]+=frac*gas_old[j];
	  Gal[p].MetalsColdGasRings[index+1]=metals_add(Gal[p].MetalsColdGasRings[index+1],gasmetal_old[j],frac);
#ifdef INDIVIDUAL_ELEMENTS
	  for (kk=0;kk<NUM_ELEMENTS;kk++)
	    Gal[p].ColdGasRings_elements[index+1][kk]+=(gasElements_old[j][kk]*frac);
#endif
	}
    }
  if (DiskRadiusModel == 0)
    Gal[p].ColdGasRadius=get_gas_disk_radius(p);

}
