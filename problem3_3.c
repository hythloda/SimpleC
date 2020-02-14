#include <stdio.h>         /*for printf*/
#include <math.h>


double newton(double x)
{
  return (x*exp(x));
}

double newtonp(double x)
{
  return(exp(x)+x*exp(x));
}

main()
{
  double delta, iterations, root, guess, guess2, midpoint;
  iterations=delta=root=delta=midpoint=0;
  delta=(newton(guess)/newtonp(guess));
  guess2=0.5;
  guess= -0.5;
  if(newton(guess)<1e-2)
    {
      for(guess=-1,root=guess-delta;fabs(newton(root))>1e-10;guess=root)
	{
	  delta=(newton(guess)/newtonp(guess));
	  root=guess-delta;
	  iterations=iterations+1;
	  printf("iterations= %lf guess= %lf\n",iterations,guess);
	}
    }
  else
    {
      while(fabs(guess-guess2)>1e-10)	{
	  midpoint=((guess + guess2)/2);
	  if(newton(guess)*newton(guess2)>0) {
	    printf("Pick different guess\n");
	  }
	  if (newton(midpoint)<0) {
	    guess=midpoint;
	  }
	  else {
	    guess2=midpoint;
	  }
	  iterations=iterations+1;
	}
    }
  printf("root= %lf\n",root);
  printf("delta= %lf\n",delta);
  printf("iterations= %lf\n",iterations);

}
