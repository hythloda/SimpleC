#include <stdio.h>       
#include <math.h>
double Integrand(double x)
{
  return (pow(x,2));
}
double func(double x)
{
  return(exp(-Integrand(x)));
}
double f(double x, double y)
{
  return(fabs(y)*(func(x)+func(x+y))/2);
}

main()
{
  double step, I, i,area , j, e, step2, I2, m, area2, total;

  j=10000;
  step=10.0/j;
  I=0.0;
  area=0.0;
  step2=-10.0/j;
  I2=0.0;
  area2=0.0;
  for(i=0.0;i<j;i=i++)
     {
       area=f(I,step)+area;
       I=I+step;
     }
printf("total area on right= %f\n",area);
  for(m=0.0;m<j;m=m++)
     {
       area2=f(I2,step2)+area2;
       I2=I2+step2;
     }
printf("total area on left= %f\n",area2);
   total=area+area2;
   e=fabs((total-1.77245385)/(total)*100);
printf("error= %f\n",e);
printf("total area= %f\n",total);
}
