#include <stdio.h>
#include <math.h>
#include "philsp.h"

int main()
{
  double xmin, xmax, ymin, ymax, theta,  radius, x[5],y[5];
  int j,u,k;
  open_plot("800x800");
  xmin=-1;
  xmax=1;
  ymin=-1;
  ymax=1;
  box_plot(xmin, xmax, ymin, ymax,1.5, 1,"x","y","","");
  for(j=1;j<10000;j++)
    {
      theta=10*(j/1000.0)*M_PI;
      radius=0.5*cos(j/100);
      for(k=0;k<5;k++)
	{
	  x[k]=radius*2.0*sin(theta*k/2.0);
	  y[k]=radius*2.0*cos(theta*k);
	}
      for(k=0;k<5;k++)
	{
	  putpoint_plot(x[k],y[k],10,1,2,2.0,0);
	}
      flush_plot();
      erase_plot();
      delay_plot(1);
    }

  return(0);
}
