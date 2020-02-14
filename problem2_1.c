#include <stdio.h>         /*for printf*/
#include <math.h>        /*for sqrt()*/

main()

{
  float x ,y;
   for(x=0.0;x<=10.0;x=x+0.25)
    {      /*no semi-colen needed*/
     y=((sqrt(x))*(cos(x)));
      /*brackat is to end the for loop*/
  printf("%4.2f %8.5f\n",x,y);
    }
  return(0);
}
