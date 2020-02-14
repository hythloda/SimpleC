#include <stdio.h>
#include <math.h>
#include <stdlib.h>
// By Amanda Martin
// equation were from Data and Error Analysis
// Does a least squares fit given points x, and y, and errors
//  prints out chi^2, int, slope, error in a, error in b

#define f(i) (slope*x[i] + intercept)
#define sqr(x) ((x)*(x))
#define s dely

int main(){
  int i, data;
  double x[6], y[6], delx[6], dely[6];
  double del, sumofx, sumofy,sumofxx,sumofyy, sumofxy, una, unb;
  double  sumofxxs, sumofxs, sumofss, sumofxys, sumofys;
  double  slope, ave, intercept, chi;
  FILE *datafile;
  datafile=fopen("steve.dat","r");
  //The name of the data file goes in the first quote.

  sumofx=sumofy=sumofxxs=sumofxs=una=unb=0;
  sumofss=sumofxys=sumofys=0;

  // !!!!**** change the value of data ****!!!!
  data=6; //This is the number of rows in the data file.
  // !!!!**** change the value of data ****!!!!

  //This is were it reads in the file.
  for(i=0;i<data;i++)
    fscanf(datafile, "%lf %lf %lf %lf\n",
	   &x[i],&y[i],&dely[i],&delx[i]);

  for(i=0;i<data;i++){
    s[i]=dely[i];
    sumofx   += x[i];
    sumofy   += y[i];
    sumofxx  += x[i]*x[i];
    sumofyy  += y[i]*y[i];
    sumofxy  += x[i]*y[i];
    sumofxxs += (x[i]*x[i])*s[i];
    sumofxs  += x[i]*s[i];
    sumofss  += s[i];
    sumofxys +=(x[i]*y[i])*s[i];
    sumofys  +=(y[i])*s[i];
  }

  del=sumofss*sumofxxs-sumofxs*sumofxs;
  ave=data*sumofxx-sumofx*sumofx;
  slope=(sumofss*sumofxys-sumofxs*sumofys)/del;
  intercept=(sumofxxs*sumofys-sumofxs*sumofxys)/del;
  una=sqrt((1/del)*(sumofxxs));
  unb=sqrt((1/del)*(sumofss));

  for(chi=0,i=0;i<data;i++)
     chi += sqr(y[i] - f(i))*s[i];
  //Seems like there is still some bug in chi.

  fclose(datafile);
  printf("chi= %lf\t slope= %e\t intercept= %e\t sigma_slope= %lf\t sigma_intercept= %lf\n",
	 chi/(data-1), slope, intercept, unb, una);

}
