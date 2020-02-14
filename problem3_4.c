#include <stdio.h>         /*for printf*/
#include <math.h>
#include <stdio.h>

double f(double x){
  return (x*x-2);
}

  double nt(double x){
  return(2*x);
}
  double ntp(double x){
  return(2);
}

int main(){
    double i, r, g, g2, m;
  i=0;
  g2=2;
  g=-1;
  /*if(nt(g)<1e-2)
      while(fabs(nt(r))>1e-10){
	  r=g-(nt(g)/ntp(g));
	  i=i+1;
	  printf("i= %6.4f g= %6.4f\n",i,g);
	  } */
  /* else */
      while(fabs(g-g2)>1e-10)	{
	  m=((g + g2)/2);
	  if(f(g)*f(g2)>0)
	    printf("Pick different g\n");
	  if (nt(m)<0)
	    g2=m;
     	  else
	    g=m;

	  i=i+1;
	}
  printf("r= %lf\n",r);

  return(0);
}
