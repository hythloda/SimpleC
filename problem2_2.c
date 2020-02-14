#include <stdio.h>         /*for printf*/
#include <math.h>        /*for cos(-1.0)*/

main()

{
  float a, b, c;
  int x, y, z;
  a=1.0;
  b=2.0;
  x=1;
  y=2;
  c=a/b;
  z=x/y;
  printf("1.0/2.0= %9.1f\n",c);
  printf("1/2= %9.1d\n",z);
return(0);
}
