#include <stdio.h>
#include <math.h>

#define NR_END 1
#define FREE_ARG char*

#define SQR(a) (a == 0.0 ? 0.0 : a*a)
#define FMAX(a,b) ((a) > (b) ? (a) : (b))
#define FMIN(a,b) ((a) < (b) ? (a) : (b))
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))


double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;

        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        if (!v) fprintf(stderr,"allocation failure in dvector()\n");
        return v-nl+NR_END;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
        long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
        double **m;

        /* allocate pointers to rows */
        m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
        if (!m) fprintf(stderr,"allocation failure 1 in matrix()\n");
        m += NR_END;
        m -= nrl;

        /* allocate rows and set pointers to them */
        m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
        if (!m[nrl]) fprintf(stderr,"allocation failure 2 in matrix()\n");
        m[nrl] += NR_END;
        m[nrl] -= ncl;

        for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

        /* return pointer to array of pointers to rows */
        return m;
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
        free((FREE_ARG) (v+nl-NR_END));
}


void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
        free((FREE_ARG) (m[nrl]+ncl-NR_END));
        free((FREE_ARG) (m+nrl-NR_END));
}


void mmid(double y[], double dydx[], int nvar, double xs, double htot, int nstep,
        double yout[], void (*derivs)(double, double[], double[]))
{
        int n,i;
        double x,swap,h2,h,*ym,*yn;

        ym=dvector(1,nvar);
        yn=dvector(1,nvar);
        h=htot/nstep;
        for (i=1;i<=nvar;i++) {
                ym[i]=y[i];
                yn[i]=y[i]+h*dydx[i];
        }
        x=xs+h;
        (*derivs)(x,yn,yout);
        h2=2.0*h;
        for (n=2;n<=nstep;n++) {
                for (i=1;i<=nvar;i++) {
                        swap=ym[i]+h2*yout[i];
                        ym[i]=yn[i];
                        yn[i]=swap;
                }
                x += h;
                (*derivs)(x,yn,yout);
        }
        for (i=1;i<=nvar;i++)
                yout[i]=0.5*(ym[i]+yn[i]+h*yout[i]);
        free_dvector(yn,1,nvar);
        free_dvector(ym,1,nvar);
}


#define MAXSTP 10000
#define TINY 1.0e-30

int kmax=0,kount;
double *xp,**yp,dxsav;

void odeint(double ystart[], int nvar, double x1, double x2, double eps, double h1,
        double hmin, int *nok, int *nbad,
        void (*derivs)(double, double [], double []),
        void (*rkqs)(double [], double [], int, double *, double, double, double [],
        double *, double *, void (*)(double, double [], double [])))
{
        int nstp,i;
        double xsav,x,hnext,hdid,h;
        double *yscal,*y,*dydx;

        yscal=dvector(1,nvar);
        y=dvector(1,nvar);
        dydx=dvector(1,nvar);
        x=x1;
        h=SIGN(h1,x2-x1);
        *nok = (*nbad) = kount = 0;
        for (i=1;i<=nvar;i++) y[i]=ystart[i];
        if (kmax > 0) xsav=x-dxsav*2.0;
        for (nstp=1;nstp<=MAXSTP;nstp++) {
                (*derivs)(x,y,dydx);
                for (i=1;i<=nvar;i++)
                        yscal[i]=fabs(y[i])+fabs(dydx[i]*h)+TINY;
                if (kmax > 0 && kount < kmax-1 && fabs(x-xsav) > fabs(dxsav)) {
                        xp[++kount]=x;
                        for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                        xsav=x;
                }
                if ((x+h-x2)*(x+h-x1) > 0.0) h=x2-x;
                (*rkqs)(y,dydx,nvar,&x,h,eps,yscal,&hdid,&hnext,derivs);
                if (hdid == h) ++(*nok); else ++(*nbad);
                if ((x-x2)*(x2-x1) >= 0.0) {
                        for (i=1;i<=nvar;i++) ystart[i]=y[i];
                        if (kmax) {
                                xp[++kount]=x;
                                for (i=1;i<=nvar;i++) yp[i][kount]=y[i];
                        }
                        free_dvector(dydx,1,nvar);
                        free_dvector(y,1,nvar);
                        free_dvector(yscal,1,nvar);
                        return;
                }
                if (fabs(hnext) <= hmin) {
                  fprintf(stderr,"Step size too small in odeint\n");
                  h = hmin;
                }
                else
                  h=hnext;
        }
        fprintf(stderr,"Too many steps in routine odeint\n");
}
#undef MAXSTP
#undef TINY


extern double **d,*x;

void pzextr(int iest, double xest, double yest[], double yz[], double dy[], int nv)
{
        int k1,j;
        double q,f2,f1,delta,*c;

        c=dvector(1,nv);
        x[iest]=xest;
        for (j=1;j<=nv;j++) dy[j]=yz[j]=yest[j];
        if (iest == 1) {
                for (j=1;j<=nv;j++) d[j][1]=yest[j];
        } else {
                for (j=1;j<=nv;j++) c[j]=yest[j];
                for (k1=1;k1<iest;k1++) {
                        delta=1.0/(x[iest-k1]-xest);
                        f1=xest*delta;
                        f2=x[iest-k1]*delta;
                        for (j=1;j<=nv;j++) {
                                q=d[j][k1];
                                d[j][k1]=dy[j];
                                delta=c[j]-q;
                                dy[j]=f1*delta;
                                c[j]=f2*delta;
                                yz[j] += dy[j];
                        }
                }
                for (j=1;j<=nv;j++) d[j][iest]=dy[j];
        }
        free_dvector(c,1,nv);
}


#define KMAXX 8
#define IMAXX (KMAXX+1)
#define SAFE1 0.25
#define SAFE2 0.7
#define REDMAX 1.0e-5
#define REDMIN 0.7
#define TINY 1.0e-30
#define SCALMX 0.1

double **d,*x;

void bsstep(double y[], double dydx[], int nv, double *xx, double htry, double eps,
        double yscal[], double *hdid, double *hnext,
        void (*derivs)(double, double [], double []))
{
        void mmid(double y[], double dydx[], int nvar, double xs, double htot,
                int nstep, double yout[], void (*derivs)(double, double[], double[]));
        void pzextr(int iest, double xest, double yest[], double yz[], double dy[],
                int nv);
        int i,iq,k,kk,km;
        static int first=1,kmax,kopt;
        static double epsold = -1.0,xnew;
        double eps1,errmax,fact,h,red,scale,work,wrkmin,xest;
        double *err,*yerr,*ysav,*yseq;
        static double a[IMAXX+1];
        static double alf[KMAXX+1][KMAXX+1];
        static int nseq[IMAXX+1]={0,2,4,6,8,10,12,14,16,18};
        int reduct,exitflag=0;

        d=dmatrix(1,nv,1,KMAXX);
        err=dvector(1,KMAXX);
        x=dvector(1,KMAXX);
        yerr=dvector(1,nv);
        ysav=dvector(1,nv);
        yseq=dvector(1,nv);
        if (eps != epsold) {
                *hnext = xnew = -1.0e29;
                eps1=SAFE1*eps;
                a[1]=nseq[1]+1;
                for (k=1;k<=KMAXX;k++) a[k+1]=a[k]+nseq[k+1];
                for (iq=2;iq<=KMAXX;iq++) {
                        for (k=1;k<iq;k++)
                                alf[k][iq]=pow(eps1,(a[k+1]-a[iq+1])/
                                        ((a[iq+1]-a[1]+1.0)*(2*k+1)));
                }
                epsold=eps;
                for (kopt=2;kopt<KMAXX;kopt++)
                        if (a[kopt+1] > a[kopt]*alf[kopt-1][kopt]) break;
                kmax=kopt;
        }
        h=htry;
        for (i=1;i<=nv;i++) ysav[i]=y[i];
        if (*xx != xnew || h != (*hnext)) {
                first=1;
                kopt=kmax;
        }
        reduct=0;
        for (;;) {
                for (k=1;k<=kmax;k++) {
                        xnew=(*xx)+h;
                        if (xnew == (*xx)) fprintf(stderr,"step size underflow in bsstep\n");
                        mmid(ysav,dydx,nv,*xx,h,nseq[k],yseq,derivs);
                        xest=SQR(h/nseq[k]);
                        pzextr(k,xest,yseq,y,yerr,nv);
                        if (k != 1) {
                                errmax=TINY;
                                for (i=1;i<=nv;i++) errmax=FMAX(errmax,fabs(yerr[i]/yscal[i]));
                                errmax /= eps;
                                km=k-1;
                                err[km]=pow(errmax/SAFE1,1.0/(2*km+1));
                        }
                        if (k != 1 && (k >= kopt-1 || first)) {
                                if (errmax < 1.0) {
                                        exitflag=1;
                                        break;
                                }
                                if (k == kmax || k == kopt+1) {
                                        red=SAFE2/err[km];
                                        break;
                                }
                                else if (k == kopt && alf[kopt-1][kopt] < err[km]) {
                                                red=1.0/err[km];
                                                break;
                                        }
                                else if (kopt == kmax && alf[km][kmax-1] < err[km]) {
                                                red=alf[km][kmax-1]*SAFE2/err[km];
                                                break;
                                        }
                                else if (alf[km][kopt] < err[km]) {
                                        red=alf[km][kopt-1]/err[km];
                                        break;
                                }
                        }
                }
                if (exitflag) break;
                red=FMIN(red,REDMIN);
                red=FMAX(red,REDMAX);
                h *= red;
                reduct=1;
        }
        *xx=xnew;
        *hdid=h;
        first=0;
        wrkmin=1.0e35;
        for (kk=1;kk<=km;kk++) {
                fact=FMAX(err[kk],SCALMX);
                work=fact*a[kk+1];
                if (work < wrkmin) {
                        scale=fact;
                        wrkmin=work;
                        kopt=kk+1;
                }
        }
        *hnext=h/scale;
        if (kopt >= k && kopt != kmax && !reduct) {
                fact=FMAX(scale/alf[kopt-1][kopt],SCALMX);
                if (a[kopt+1]*fact <= wrkmin) {
                        *hnext=h/fact;
                        kopt++;
                }
        }
        free_dvector(yseq,1,nv);
        free_dvector(ysav,1,nv);
        free_dvector(yerr,1,nv);
        free_dvector(x,1,KMAXX);
        free_dvector(err,1,KMAXX);
        free_dmatrix(d,1,nv,1,KMAXX);
}
#undef KMAXX
#undef IMAXX
#undef SAFE1
#undef SAFE2
#undef REDMAX
#undef REDMIN
#undef TINY
#undef SCALMX
#undef NRANSI
