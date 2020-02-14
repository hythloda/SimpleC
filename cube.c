#include <stdio.h>
#include <math.h>

int main()
{
  float input, cube;
  int j;
  input=4.5;
  cube=1.0;
  for(j=1;j<=3;j=j+1){
    cube=cube*input;
}
  /* we've done the cube*/
  printf("The cube of %f is %f!\n",input,cube);
  return;
}
