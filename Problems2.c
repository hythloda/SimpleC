#include <stdio.h>
#include <math.h>
#include <stdlib.h>


/*Decay solves for time to get a lot of fe

Amanda Martin October 6, 2002*/

void func( double (*new_i)(double i, double j, double k),
	   double (*new_j)(double i, double j, double k),
	   double (*new_k)(double i, double j, double k),
	   double *i, double *j, double *k, double p)
{
  while(p<50)
    {
      int o,m,l;
      o=(int)i;
      m=(int)j;
      l=(int)k;
      *i=new_i(*i,*j,*k);
      *j=new_j(*i,*j,*k);
      *k=new_k(*i,*j,*k);
      p=p+1;
      if(l%1==0, m%1==0, o%1==0)
	{
      printf("i= %d\t j= %d\t k=%d\t p=%d\n",*i,*j,*k, p);
	}
    }
}

double new_i(double i, double j, double k)
{
  i=sqrt(fabs(k*k-j*j));
  return(i);
}

 double new_j(double i,double j, double k)
{
  j=sqrt(fabs(k*k-i*i));
  return(j);
}

double new_k(double i, double j, double k)
{
  k=sqrt(fabs(i*i+j*j));
     return(k);
	}
  int main()
    {
      double  i, j, k;
      int t, p;
      t=1;
      p=0.0;
      i=1.0;
      j=1.0;
      k=1.0;
      func(&new_i, &new_j,&new_k,&i,&j,&k,p);
    }
