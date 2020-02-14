#include <stdio.h>         
#include <stdlib.h>
#include <math.h>

#define theta .0001
#define b -.0569
#define r 6000.
/*Bisection by Amanda Morrow September 2002 modifies for astronomy 302*/

double f(double x)
{
  return  (x*x*x - (2.*x*x*b)/(4.*theta)+((b*r)/4.-r*r/(2.*theta*theta))*x+(r*r*r)/(4.*theta));
}
int main(){
  double low, high, midpoint, root, iterations;
  low=-1000;
  high=10000;
  iterations=0;
  while(high-low>1e-10){
      midpoint=((low + high)/2);
      if(f(low)*f(high)>0) {
	printf("Pick different roots\n");
      }
      if (f(midpoint)<0){
	low=midpoint;
      }
      else {
	high=midpoint;
      }
      iterations=iterations+1;
      }

root= (low+high)/2;
printf("the root is= %lf\n",root);
}
