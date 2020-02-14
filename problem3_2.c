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
  double delta1, delta2, delta3, iterations1, iterations2, iterations3, root1, root2, root3, guess1, guess2, guess3;
  delta1=(newton(guess1)/newtonp(guess1));
  delta2=(newton(guess2)/newtonp(guess2));
  delta3=(newton(guess3)/newtonp(guess3));
  iterations1=0;
  iterations2=0;
  iterations3=0;
  guess1=2;
  guess2=7;
  guess3=16;
  root1=guess1-delta1;
  root2=guess2-delta2;
  root3=guess3-delta3;
  for(guess1=2;fabs(root1-guess1)>1e-5;guess1=guess1+delta1)
    {
      delta1=(newton(guess1)/newtonp(guess1));
      iterations1=iterations1+1;
    }
  for(guess2=2;fabs(root2-guess2)>1e-5;guess2=guess2+delta2)
    {
      delta2=(newton(guess2)/newtonp(guess2));
      iterations2=iterations2+1;
    }
  for(guess3=2;fabs(root3-guess3)>1e-5;guess3=guess3+delta3)
    {
      delta3=(newton(guess3)/newtonp(guess3));
      iterations3=iterations3+1;
    }
  root1=guess1+delta1;
  root2=guess2+delta2;
  root3=guess3+delta3;
  printf("root1= %lf root2= %lf root3= %lf\n",root1,root2,root3);
  printf("delta1= %lf delta2= %lf delta3= %lf\n", delta1, delta2, delta3);
}
