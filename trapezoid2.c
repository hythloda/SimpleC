#include <stdio.h>       
#include <math.h>
double Integrand(double x)
{
  return ((x*x*x)/(exp(x)-1));
}
double f(double x, double y)
{
  return(y*(Integrand(x)+Integrand(x+y))/2);
}

main()    /*takes a while but it does end*/
{
  double step, I, i,area ,j;

  j=10000;
  step=100/j;
  I=1e-6;
  area=0.0;
  for(i=0.0;i<j;i=i++)
     {
       area=f(I,step)+area;
       I=I+step;

     }
printf("step= %f area = %f i= %f\n",step,area,i);
}
