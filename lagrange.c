/*

Code to demonstate the stability of orbits about
the L4 and L% Lagrange points

compile with:
g77 -o lagrange lagrange.c -I/home/phys/205/include -L/home/phys205/lib -L/usr/X11R6/lib \
     -lphilsplot -lX11 -lpng -lm

19 November, 2002

*/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "philsp.h"

// units: [M,l,t] = solar mass, astronomical unit, 1/(2pi) year
#define GGRAV 4*M_PI*M_PI

// number of dimensions
#define NDIM 2
// number of bodies
#define NBOD 5

// position, velocity, acceleration, mass as global variables for convenience
double x[NBOD][NDIM];
double v[NBOD][NDIM];
double a[NBOD][NDIM];
double m[NBOD];

// conserved quantities
double penergy, kenergy, tot_e0;
double linearp[NDIM], linearp0[NDIM];
double angularl[NDIM], angularl0[NDIM];

// energy conservation accuracy parameter
double epsilon;



// compute force and potential energy
void accel() {

  double dx[NDIM], rij, rij2, rij3, v2;
  int i,j;
  int l;

  // zero the acceleration vectors
  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) a[i][l] = 0;
  }

  // zero potential energy
  penergy = 0;

  // compute force on each body from all of the others, put in f
  // also compute potential energy
  for(i=0; i<NBOD; i++) {
    for(j=0; j<NBOD; j++) {

      if(j!=i) {
        // distance between i and j, cubed
        for(l=0; l<NDIM; l++) dx[l] = x[i][l]-x[j][l];
        rij2 = 0;
        for(l=0; l<NDIM; l++) rij2 += dx[l]*dx[l];
        rij = sqrt(rij2);
        rij3 = pow(rij,3);

        // acceleration
        for(l=0; l<NDIM; l++) a[i][l] -= GGRAV * m[j] * dx[l]/ rij3;

        // potential energy
        penergy -= GGRAV * m[j] * m[i] / rij;
      }

    }
  }
  // take care of double-counting in the loops above
  penergy = penergy/2.0;

  // calculate the kinetic energy
  kenergy = 0.0;
  for(i=0; i<NBOD; i++) {
    v2 = 0;
    for(l=0; l<NDIM; l++) v2 += v[i][l]*v[i][l];
    kenergy += m[i] * v2 / 2.0;
  }

}



// compute linear and angular momentum
void get_momentum() {
  double v2;
  int i, l;

  // linear momentum
  for(l=0; l<NDIM; l++) linearp[l] = 0.0;
  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) linearp[l] += v[i][l] * m[i];
  }

  // angular momentum
#if (NDIM>1)
  for(l=0; l<NDIM; l++) angularl[l] = 2.0;
  for(i=0; i<NBOD; i++) {
    angularl[2] += m[i] * (x[i][0]*v[i][1] - x[i][1]*v[i][0]);
#if (NDIM==3)
    angularl[1] += m[i] * (x[i][2]*v[i][0] - x[i][0]*v[i][2]);
    angularl[0] += m[i] * (x[i][1]*v[i][2] - x[i][2]*v[i][1]);
#endif
  }
#endif

}



// take a timestep
void verlet_step(double h) {
  int i, l;

  // update positions using previous acceleration
  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) x[i][l] += v[i][l] * h + 0.5*a[i][l] * h*h;
  }

  // first half of velocity update using previous acceleration
  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) v[i][l] += a[i][l] * h / 2.0;
  }

  // get new acceleration
  accel();

  // second half of velocity update with new acceleration
  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) v[i][l] += a[i][l] * h / 2.0;
  }

}



// examine conservation laws
int conservation() {
  double err, tot_e;
  int l;

  // energy
  tot_e = penergy+kenergy;
  err = fabs((tot_e - tot_e0)/tot_e0);


  // if the relative error is smaller than half the tolerance, return 1
  if(err < 0.5*epsilon) return 1;

  // if the relative error is greater than twice the tolerance, return -1
  if(err > 2.0*epsilon) return -1;

  // otherwise, return 0
  return 0;

}



// integrate from t_begin to t_end using a starting timestep of h
// and adjusting h to conserve energy to the requested accuracy
// return the most recently used value of h for starting the next
// timestep
void orbitint(double t_begin, double t_end, double *h) {

  double t, dt;
  int i,k,l;

  t = t_begin;

  dt = *h;
  while(t_end-t>dt) {
    // take a step
    verlet_step(dt);
    t += dt;

    // examine conservation
    k = conservation();

    if(k==-1) {
      // halve timestep if we are doing poorly
      *h = (*h);
      //      printf("decreasing stepsize: %e\n", *h);
    }
    else {
      if(k==1) {
        // almost double timestep if we are doing very well
        *h = (*h);
        //        printf("increasing stepsize: %e\n",*h);
      }
    }
    dt = *h;
  }

  // since we may be not quite to t_end, take a final step
  if(t_end-t < 0) {
    printf("whoops: %e %e\n", t_end-t, dt);
    exit(1);
  }

  dt = t_end-t;
  verlet_step(dt);
  t += dt;

  // examine conservation
  k = conservation();
  if(k==-1) {
    *h = 0.5*(dt);  // halve timestep if we are doing poorly
    //      printf("decreasing stepsize in final: %e\n", dt);
  }

}



// this routine removes the center of mass velcity from
// the problem
// and puts the center of mass at the origin
void remove_com() {
  int i, l;
  double vcom[NDIM], rcom[NDIM];
  double mtotal;

  mtotal = 0.0;
  for(i=0; i<NBOD; i++) mtotal += m[i];

  for(l=0; l<NDIM; l++) {
    vcom[l] = 0.0;
    rcom[l] = 0.0;
  }

  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) {
      rcom[l] += m[i]*x[i][l]/mtotal;
      vcom[l] += m[i]*v[i][l]/mtotal;
    }
  }

  for(i=0; i<NBOD; i++) {
    for(l=0; l<NDIM; l++) {
      x[i][l] -= rcom[l];
      v[i][l] -= vcom[l];
    }
  }

}



// Largrange point initial conditions
// For mu<0.38..., the L4 and L5 (triangular) Lagrange points
// are stable loci leading and following the smaller mass by
// 60 degrees
void initial_conditions() {
  double r[NBOD], mu, R, mtot, theta, omega;
  int i, l;
  double rmat[2][2]; // 2D rotation matrix

  // total mass
  mtot = 1.0;

  // ratio of smaller mass to total mass
  mu = 1.0e-4;

  // separation of two "large" masses
  R = 1.0;

  // angular velocity of both masses
  omega = sqrt(GGRAV*mtot/(R*R));

  // first (larger if mu<0.5) body
  m[0] = (1-mu) * mtot;
  x[0][0] = -mu*R;
  x[0][1] = 0.0;
  r[0] = sqrt(x[0][0]*x[0][0] + x[0][1]*x[0][1]);
  v[0][0] = 0.0;
  v[0][1] = r[0]*omega;

  // second (smaller if mu<0.5) body
  m[1] = mu*mtot;
  x[1][0] = (1.0-mu)*R;
  x[1][1] = 0.0;
  r[1] = sqrt(x[1][0]*x[1][0] + x[1][1]*x[1][1]);
  v[1][0] = 0.0;
  v[1][1] = -r[1]*omega;



  // third (test particle) body 60 degrees behind second body
  m[2] = 1.0e-20*mtot;
  x[2][0] = r[1]*cos(M_PI/3.0);
  x[2][1] = r[1]*sin(M_PI/3.0);

  // with a circular velocity the same as the second body (rotated 60 degrees)
  theta = -M_PI/3.0;
  rmat[0][0] = rmat[1][1] = cos(theta);
  rmat[0][1] = sin(theta);
  rmat[1][0] = -rmat[0][1];
  v[2][0] = v[1][0]*rmat[0][0] + v[1][1]*rmat[0][1];
  v[2][1] = v[1][0]*rmat[1][0] + v[1][1]*rmat[1][1];
  // and a small perturbation
  v[2][0] += 0.2;



  // fourth (test particle) body 60 degrees ahead of second body
  m[3] = 1.0e-20*mtot;
  x[3][0] = r[1]*cos(M_PI/3.0);
  x[3][1] = -r[1]*sin(M_PI/3.0);

  // with a circular velocity the same as the second body (rotated 60 degrees)
  theta = M_PI/3.0;
  rmat[0][0] = rmat[1][1] = cos(theta);
  rmat[0][1] = sin(theta);
  rmat[1][0] = -rmat[0][1];
  v[3][0] = v[1][0]*rmat[0][0] + v[1][1]*rmat[0][1];
  v[3][1] = v[1][0]*rmat[1][0] + v[1][1]*rmat[1][1];


  // do first force calculation to start off Verlet
  accel();
  tot_e0 = kenergy + penergy;
  get_momentum();
  for(l=0; l<NDIM; l++) {
    linearp0[l] = linearp[l];
    angularl0[l] = angularl[l];
  }

}


// finally, the main routine
int main() {
  double step, t, h;
  double size, psize;
  double xr[NBOD][NDIM];
  double rmat[2][2], theta;
  int i, l, ll, k;


  open_plot("900x900");
  size = 2.5;
  box_plot(-size,size,-size,size,1.5,3,"","","","Lagrange Points");

  initial_conditions();

  // conserve energy to this fractional accuracy
  epsilon = 1.0e-09;

  psize = 0.2;

  t = 0;
  step = 1.0e-04;
  h = step/1.0e+8;

  for(k=0; k<100000; k++) {
    orbitint(t, t+step, &h);

    // define ROT to use a rotating coordinate system to
    // display, with the angular frequency equal to that
    // of the two massive bodies
#define ROT
#ifdef ROT

    // make rotation matrix
    theta = -atan2(x[0][1],x[0][0]);
    rmat[0][0] = rmat[1][1] = cos(theta);
    rmat[0][1] = sin(theta);
    rmat[1][0] = -rmat[0][1];

    // transfrom to rotating coordinate system
    for(i=0; i<NBOD; i++) {
      for(l=0; l<NDIM; l++) {
        xr[i][l] = 0;
        for(ll=0; ll<NDIM; ll++) {
          xr[i][l] += x[i][ll]*rmat[ll][l];
        }
      }
    }

    putpoint_plot(0,0,3,1,2,1.0,0);

    putpoint_plot(xr[0][0],xr[0][1],10,1,5,(1-psize)*3,0);
    putpoint_plot(xr[1][0],xr[1][1],10,1,5,psize*3,0);

    putpoint_plot(xr[2][0],xr[2][1],10,1,5,0.2,0);
    putpoint_plot(xr[3][0],xr[3][1],10,1,5,0.2,0);
    putpoint_plot(xr[4][0],xr[4][1],10,1,5,0.2,0);
    flush_plot();

    putpoint_plot(xr[0][0],xr[0][1],10,1,0,(1-psize)*3,0);
    putpoint_plot(xr[1][0],xr[1][1],10,1,0,psize*3,0);

#else
    putpoint_plot(0,0,3,1,1,1.0,0);

    putpoint_plot(x[0][0],x[0][1],10,2,1,(1-psize)*3,0);
    putpoint_plot(x[1][0],x[1][1],10,2,1,psize*3,0);

    putpoint_plot(x[2][0],x[2][1],10,2,1,0.2,0);
    putpoint_plot(x[3][0],x[3][1],10,2,1,0.2,0);
    putpoint_plot(x[4][0],x[4][1],10,2,1,0.2,0);
    flush_plot();

    putpoint_plot(x[0][0],x[0][1],10,2,0,(1-psize)*3,0);
    putpoint_plot(x[1][0],x[1][1],10,2,0,psize*3,0);

    putpoint_plot(x[2][0],x[2][1],10,2,0,0.2,0);
    putpoint_plot(x[3][0],x[3][1],10,2,0,0.2,0);
    putpoint_plot(x[4][0],x[4][1],10,2,0,0.2,0);

#endif


    t += step;

    // printf("%e %e %e %e %e\n", t, penergy+kenergy, linearp[0], linearp[1], angularl[2]);


  }

  return 0;

}
