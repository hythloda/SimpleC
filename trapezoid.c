#include <stdio.h>       
#include <math.h>
double Integrand(double x)
{
  return (x*exp(-x));
}
double f(double x, double y)
{
  return(y*(Integrand(x)+Integrand(x+y))/2);
}

main()
{
  double step, I, i,area ,j,e;

  j=10000;
  step=1.0/j;
  I=0.0;
  area=0.0;
  for(i=0.0;i<j;i=i++)
     {
       area=f(I,step)+area;
       I=I+step;

     }
  e=fabs((area-.26424112)/(area)*100);
printf("step= %f area = %f i= %f\n",step,area,i);
printf("error= %f\n",e);
}
