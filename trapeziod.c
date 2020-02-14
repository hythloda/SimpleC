#include <stdio.h> 
#include <math.h>
double f(double x)
{
return (x*exp(-x));
}
int main ()
{
  double a,b,c,d,g,h,k;
  a=0.0;/*start=a*/
  b=1.0;/*end=b*/
  c=1.0;/*step=c*/
  d=1.0;/*scaling=d*/
  while(c>=(1.0/16.0))
    {
      h=((b-a)/2)*(f(a)+f(b));/*one time through*/
      k=((b-a)/2)*(f(a)+f(c)+f(b-c)+f(b-(2.0)*c)+f(b-(3.0)*c)+f(b-(5.0)*c)+f(b-(5.0)*c)+f(b-(6.0)*c)+f(b-(7.0)*c)+f(b-(8.0)*c)+f(b-(9.0)*c)+f(b-(10.0)*c)+f(b-(11.0)*c)+f(b-(12.0)*c)+f(b-(13.0)*c)+f(b-(14.0)*c)+f(b-(15.0)*c)+f(b-(16.0)*c));
      g=h+k;/*area=g*/
      c=c/(2.0);
    }
printf("k= %f\n",k);
}
