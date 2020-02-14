#include <stdio.h>         /*for printf*/
#include <math.h>


double newton(double x)
{
  return (x*x*x-25*x*x+165*x-275);
}

double newtonp(double x)
{
  return(3*x*x-50*x+165);
}

main()
{
  double delta, iterations, root, guess;
 /*three diffrent guesses should be near 2,7, and 15*/
  iterations=0;
  delta=(newton(guess)/newtonp(guess));
 for(guess=15,root=guess-delta;fabs(newton(root))>1e-10;guess=root)
    {
      delta=(newton(guess)/newtonp(guess));
      root=guess-delta;
      iterations=iterations+1;
      printf("iterations= %4.0f guess= %6.4f\n",iterations,guess);
    }
  printf("root= %lf\n",root);
  printf("delta= %lf\n",delta);
  printf("iterations= %4.0f\n",iterations);
}
