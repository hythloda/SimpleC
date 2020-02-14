#include <stdio.h>         /*for printf*/
#include <math.h>        /*for pow and acos(-1.0)*/
#include <stdlib.h>    /*for drand48*/

int main()

{

  int ndarts,i,cnt;
  float x,y,r2,pi;
  srand48(62976297);   /*for the seed*/
  ndarts=1000000;
  cnt=0;
  for(i=0;i<ndarts;i=i+1)
    {      /*no semi-colen needed*/
      x=drand48();
      y=drand48();
      r2=pow(x,2)+pow(y,2);
      if(r2<=1.0)cnt=cnt+1;
    }   /*brackat is to end the for loop*/
  pi=4.0*(float)cnt/(float)ndarts;
  printf("pi is %f\n",pi);
  printf("really %f\n",acos(-1.0));
  return(0);
}
