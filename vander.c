#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "philsp.h"

/*Amanda Morrow
 November 11, 2002
 This takes a differential equation (van der pol) and shows a graph of the velocity vs the position.  Go down to main to change the initial conditions.  It is taken here that x is the velocity and y is the position.  Consider this like a pendulum where A is amplitute, w is hte frequency of ossiclation, i is the total time passed, and lambda is the dambing force.*/

void func( double (*first)(double x, double  y, double dt, double i),
	     double (*second)(double x, double y, double dt),
	     double *x, double *y, double dt, double *i)
{
  while(*i<10000)
    {
      *x=*x+first(*x, *y, dt, *i);
      *y=*y+second(*x, *y, dt);
      *i=*i+dt;
      putpoint_plot(*x,*y,2,1,3,2.0,0);
      flush_plot();
      delay_plot(1);

    }
}
double second( double x, double y,double dt)
{
  y=x*dt;
  return(y);
}
double first(double x, double  y, double dt, double i)
{
  double w, lambda, A;
  w=.1;
  A=.1;
  lambda=1;
  x=(A*cos(w*i)-y-lambda*(1-y*y)*x)*dt;
  return(x);
}

int main()
{

 double xmin, xmax, ymin, ymax, dt, x, y, i, w, lambda, A;
 int j,u,k;
 open_plot("800x800");
 xmin=-10;
 xmax=10;
 ymin=-10;
 ymax=10;
 dt=0.1;
 x=.20;
 y=2.0;
 i=0.0;
 box_plot(ymin, ymax, xmin, xmax,1.5, 1,"Position","velocity","","");
 func( &first, &second, &x, &y, dt, &i);
 return(0);
}
