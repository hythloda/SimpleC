#include <stdio.h>
#include <math.h>

double dif(double T)
{
  return(1/((.002)*(T-24)));
}
double sol(double x, double y)
{
  return(y*(dif(x)+dif(x+y))/2);
}

main()
{
  double step, I,  T , e;
  step=0.0001;
  I=0.0;
  T=100.0;
  while(I<900)
     {
       I=sol(T,step)+I;
       T=T-step;
     }
  e=fabs(((T-36.5627)/(T))*100);
  printf(" step= %f\n Temperature after 15 minutes= %f\n I= %f\n error= %f\n",step,T,I,e);
}
