#include <stdio.h>
#include <math.h>

/*given two values of the dependant and independant this solves for the constant.
  Amanda Morrow October 7,2002*/

double dif(double T)
{
  return(1/(T-24));
}

double sol(double x, double y)
{
  return(y*(dif(x)+dif(x-y))/2);
}
main()
{
  double step, I, T,K,e;
  int i;
  step=0.0001;
  I=0.0;
  T=100.0;
  i=0;
  while(T>33.2)
     {
       I=sol(T,step)+I;
       T=T-step;
       i=i+1;
      }

       K=I/(900);
  e=fabs(((K-.0023461442847014)/(K))*100);
       printf(" the value for K is= %lf I= %4.0f the error is %f i= %d\n",K,I,e,i);
}
