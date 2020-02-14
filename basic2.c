#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "philsp.h"

// We will use units of
// Astronomical Units for length
// Solar Masses for mass
// 1/(2 Pi) Earth Year for time
#define GGRAV 4*M_PI*M_PI

// NBODMX is the maximum number of bodies (masses) for which
// to provide space
#define NBODMX 403

// NDIM is the number of spatial dimensions to use
#define NDIM 2

// MAP is a function to map x[i][l] and v[i][l] to the
// single, one-dimensional, one-offset (y[1..N]) array
// used by odeint, below
// k=0 corresponds to position, k=1 to velocity
// l=0 is X, l=1 is Y, and l=2 is Z direction (l=2 only if NDIM is 3)
// i is the body number
// Given a set (i,l,k), MAP(i,l,k) produces a unique index into the array
#define MAP(i,l,k) (1+l*(NDIM-1)+k*NDIM+i*NDIM*2)

// these are mnemonics to use with MAP
// so, for example, the Y velocity of body i is foo[MAP(i,Y,V)]
#define X 0
#define Y 1
#define Z 2
#define R 0
#define V 1

// number of masses
int nbodies;
// position and velocity
double x[NBODMX][NDIM], v[NBODMX][NDIM];
// mass
double m[NBODMX];
// semimajor axis, eccentricity of orbits
double semi[NBODMX], eccen[NBODMX];
// energy variables
double energy, penergy, kenergy;
// total angular momentum vector
double l_tot[3];
// function prototypes to use odeint routine below

void bsstep(double y[], double dydx[], int nv, double *xx, double htry,
            double eps, double yscal[], double *hdid, double *hnext,
            void (*derivs)(double, double [], double []));
// In the following code, loop variables i and j are used for body number
// l is used for direction
// compute derivatives
// odeint gives us the independent variable (here called t),
// and the values of the dependent variables (here called f)
// it expects us to compute the derivatives of these variables
// and return them in the array df
// f and df are one-offset, one-dimensional vectors
// we use the MAP definition above to translate between this and
// a more meaningful notation

void derivs(double t, double *f, double *df)
{
  int i, j, l;
  double dr[NDIM];
  double rij, rij2, rij3;

  // zero the derivatives
  for(i=0; i<nbodies; i++)
    {
    for(l=0; l<NDIM; l++)
      {
	df[MAP(i,l,R)] = 0;   // position derivative w/r to time
	df[MAP(i,l,V)] = 0;   // velocity derivative w/r to time
      }
    }

  // zero potential energy
  penergy = 0;

  // Loop over all bodies to compute the derivatives
  for(i=0; i<nbodies; i++) {

    // compute velocity derivative
    // sum the acceleration due to all other bodies j
    for(j=0; j<nbodies; j++) {

      // self-acceleration not allowed!
      if(i!=j) {
        // get displacement vector
        for(l=0; l<NDIM; l++) dr[l] = f[MAP(i,l,0)] - f[MAP(j,l,0)];
        // get distance
        rij2 = 0;
        for(l=0; l<NDIM; l++) rij2 += dr[l]*dr[l];
        rij = sqrt(rij2);
        rij3 = rij*rij*rij;
        // the acceleration of body i due to body j
        for(l=0; l<NDIM; l++) df[MAP(i,l,V)] -= GGRAV * m[j] *dr[l] / rij3;
        // while we have the data handy, compute the potential energy
        penergy += GGRAV*m[i]*m[j]/rij;
      }
    }

    // compute position derivative (aka velocity!)
    for(l=0; l<NDIM; l++) df[MAP(i,l,R)] = f[MAP(i,l,V)];
  }
  // the algorithm above double-counts (i.e. adds
  // the constribution to the energy twice for each pair of bodies
  penergy = penergy/2.0;

}
// integrate the orbits
// from time t1 to time t2 starting with timestep hguess,
// keeping the solution relative accuracy to within tol
void orbit_int(double t1, double t2, double hguess, double tol) {
  static double y[NBODMX*NDIM*2];
  int i, l;
  int nok, nbad;
  double hmin;

  // odeint requires the minimum timestep allowed (can be zero)
  hmin = 1.0e-06*hguess;

  // load the initial conditions into the 1-based array for odeint
  for(i=0; i<nbodies; i++) {
    for(l=0; l<NDIM; l++) {
      y[MAP(i,l,R)] = x[i][l];
      y[MAP(i,l,V)] = v[i][l];
    }
  }

  // integrate the nbodies*NDIM*2 equations from t1 to t2
  odeint(y, nbodies*NDIM*2, t1, t2, tol, hguess, hmin, &nok, &nbad, derivs, bsstep);

  // now put the result back into our x and v arrays
  for(i=0; i<nbodies; i++) {
    for(l=0; l<NDIM; l++) {
      x[i][l] = y[MAP(i,l,R)];
      v[i][l] = y[MAP(i,l,V)];
    }
  }
}

// removes the center of mass velcity from the system
// and put the center of mass at the origin
void remove_com() {
  int i, l;
  double vcom[NDIM], rcom[NDIM];
  double mtotal;

  // get total mass
  mtotal = 0.0;
  for(i=0; i<nbodies; i++) mtotal += m[i];

  // zero in preparation for summing
  for(l=0; l<NDIM; l++) {
    vcom[l] = 0.0;
    rcom[l] = 0.0;
  }

  // find center of mass position and velocity
  for(i=0; i<nbodies; i++) {
    for(l=0; l<NDIM; l++) {
      rcom[l] += m[i]*x[i][l]/mtotal;
      vcom[l] += m[i]*v[i][l]/mtotal;
    }
  }
  // subtract from the bodies
  for(i=0; i<nbodies; i++) {
    for(l=0; l<NDIM; l++) {
      x[i][l] -= rcom[l];
      v[i][l] -= vcom[l];
    }
  }
}
void initial_conditions() {
  int i;
  double  r, dr,K, rmin, rmax, body,mtot;
  nbodies = 2;
  body=nbodies;
  K=.2;
  mtot=1;
  m[0] = 1/body;
  x[0][0] = 1.0;
  x[0][1] = 0.0;
  v[0][0] = 0.0;
  v[0][1] = sqrt(K*GGRAV*(mtot)/x[1][0]);

  for(i=1; i<nbodies; i++) {
    m[i] = 1/body;
    x[i][0] = cos(2*M_PI/i);
    x[i][1] = sin(2*M_PI/i);
    v[i][0] = -sqrt(K*GGRAV*(mtot)/x[i][1])*sin(2*M_PI/i);
    v[i][1] = sqrt(K*GGRAV*(mtot)/x[i][0])*cos(2*M_PI/i);
    remove_com();
  }
}

void conservation() {
  int i, l;
  double v2;

  // we have potential energy from derivs, above
  // get kinetic energy here
  kenergy = 0.0;
  for(i=0; i<nbodies; i++) {
    v2 = 0.0;
    for(l=0; l<NDIM; l++) v2 += v[i][l]*v[i][l];
    kenergy += m[i]*v2/2.0;
  }
  // total energy
  energy = kenergy + penergy;
  // get angular momentum
  // for NDIM of 1, there is no angular momentum
  for(l=0; l<3; l++) l_tot[l] = 0.0;
#if (NDIM>1)
  // for NDIM of 2, the angular momentum vector points out of the plane
  //    so we need only one component
  for(i=0; i<nbodies; i++) {
    l_tot[2] += m[i] * (x[i][0]*v[i][1] - x[i][1]*v[i][0]);
#if (NDIM==3)
    // for NDIM of 3, we need all three components of the angular momentum
    l_tot[1] += m[i] * (x[i][2]*v[i][0] - x[i][0]*v[i][2]);
    l_tot[0] += m[i] * (x[i][1]*v[i][2] - x[i][2]*v[i][1]);
#endif
  }
#endif

  // for the two-dimensional version
  // printf("conservation: %12.5e  %12.5e\n", energy, l_tot[2]);
}

int main() {
  int i, k, n, nplot;
  double t, step, hguess, tol;
  double size, psize;
  open_plot("800x800");  // set up the plotting window
  psize = 0.5;           // the size of the points to plot
  nplot = 1;             // plot every nplot calls to odeint
  n = 0;
  size = 6.9;            // the limits of the box
  box_plot(-size,size,-size,size,1.5,1,"","","","");
  flush_plot();
  // set up the problem
  t = 0;
  initial_conditions();
  step = 1.0e-2;           // each call to odeint advances by this time
  hguess = step/1.0e+3;    // start with this timestep inside odeint
  tol = 1.0e-03;           // keep solution relative accuracy to this

  // loop over time
  for(k=0; k<100000; k++) {
    // integrate from t to t+step
    orbit_int(t, t+step, hguess, tol);
    t += step;
    n++;
    if(n==nplot) {
      n = 0;
      // if(k%100==0){
      //	sleep(1);
      //}
      for(i=0; i<nbodies; i++) {
        putpoint_plot(x[i][0],x[i][1],10,1,5,psize,0);
      }
      flush_plot();
      for(i=0; i<nbodies; i++) {
        putpoint_plot(x[i][0],x[i][1],10,1,0,psize,0);
      }
       conservation();
    }
  }
  return 0;
}
