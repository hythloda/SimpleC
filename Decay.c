#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/*Decay solves for time to get a lot of fe

Amanda Morrow October 6, 2002*/

double new_ni(double ni, double t)
{
  return(ni*(1-(t/8.8)));
}
double new_co(double co, double ni, double t)
{
  return(co*(1-(t/111))+(ni*(t/8.8)));
}
double new_fe(double co, double ni)
{
  return(1-ni-co);
}

int main()
{
  double  ni, co, fe, dt;
  int t;
  t=0;
  ni=1.0;
  co=0.0;
  fe=0.0;
   while(fe<.99)
    {
      ni=fabs(new_ni(ni,t));
      co=new_co(co,ni,t);
      fe=new_fe(co,ni);
      t=t+1;

    }
  printf("ni= %f\t co= %f\t fe= %f\t the time elapsed is %d days\n",ni, co, fe, t);
}
