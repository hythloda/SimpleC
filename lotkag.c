#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "philsp.h"

/*Amanda Morrow
 November 11, 2002
This takes a differential equation (lotka-Voltera) and shows a graph.the program lotka.c is the one that will output the numbers for an initial rabit population.   dt, R, F, i declared in main.  Change these values when the initial conditions need to be changed.   */


void func( double (*first)(double R, double  F, double dt),
	     double (*second)(double R, double F, double dt),
	     double *R, double *F, double dt, double *i)
{
  while(*i<10000)
    {
      *R=*R+first(*R, *F, dt);
      *F=*F+second(*R, *F, dt);
      *i=*i+dt;
      putpoint_plot(*R,*F,2,1,8,2.0,0);
      flush_plot();
      delay_plot(1);
    }
}

double first(double R, double  F, double dt)
{
  double a,b;
  a=0.2;
  b=.01;
  R=(a*R-b*R*F)*dt;
  return(R);
}
double second( double R, double F,double dt)
{
  double c,d;
  c=.001;
  d=.1;
  F=(c*F*R-d*F)*dt;
  return(F);
}

int main()
{
 double Rmin, Rmax, Fmin, Fmax, dt, R, F, i;
 int j,u,k;
 open_plot("800x800");
 Rmin=0;
 Rmax=300;
 Fmin=0;
 Fmax=50;
 dt=0.1;
 R=200.0;
 F=20.0;
 i=0.0;
 box_plot(Rmin, Rmax, Fmin, Fmax,1.5, 1,"Rabits","Foxes","","");
 func(&first,&second, &R, &F, dt, &i);
 return(0);
}
