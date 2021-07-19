// fdist.cpp : Defines the entry point for the console application.
//

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nmath.h"
#include "dpq.h"

#include "stat.h"

void nrerror(char* msg) {
  printf("ERR: %s\n", msg);
}


double betacf(double a, double b, double x) {

  /* Used by betai: evaluates continued fraction for incomplete beta
     function by modified Lentz's method */

 int m, m2;
 double aa,c,d,del,h,qab,qam,qap;

 qab=a+b;
 qap= a+1.0;
 qam= a-1.0;
 c = 1.0;
 d = 1.0 - qab*x/qap;
 if (fabs(d)< FPMIN) d=FPMIN;
 d = 1.0 /d;
 h = d;
 for (m=1; m<MAXIT; m++) {
  m2 = 2*m;
  aa = m*(b-m)*x / ((qam+m2)*(a+m2));
  d=1.0+aa*d;
  if (fabs(d) < FPMIN) d = FPMIN;
  c=1.0+aa/c;
  if (fabs(c) <FPMIN) c = FPMIN;
  d = 1.0/d;
  h *= d*c;
  aa = -(a+m)*(qab+m)*x / ((a+m2)*(qap+m2));
  d = 1.0 + aa*d;
  if (fabs(d) < FPMIN) d=FPMIN;
  c = 1.0+ aa/c;
  if (fabs(c) <FPMIN) c=FPMIN;
  d = 1.0/d;
  del = d*c;
  h *= del;
  if (fabs(del-1.0) <EPS) break;
 }
 if (m > MAXIT) nrerror("a or b too big, or MAXIT too small in betacf");
 return h;
}
double gammln(double xx)
     /* returns the value ln[gamma(xx)] for xx > 0 */
{
  /*
   Internal arithmetic will be done in double precision, a nicety that you
   can omit if five-figure accuracy is good enough. */

 double x,y, tmp, ser;
 static double cof[6]= {76.18009172947146,-86.50532032941677,
  24.01409824083091, -1.231739572450155,0.1208650973866179e-2, -0.5395239384953e-5};
 int j;

 y=x=xx;
 tmp = x+5.5;
 tmp -= (x+0.5)*log(tmp);
 ser = 1.000000000190015;
 for (j=0; j<= 5; j++) ser += cof[j]/++y;
 return -tmp+log(2.5066282746310005*ser/x);
}


double betai(double a, double b, double x)
{
  double bt;
  if (x < 0.0 || x > 1.0) nrerror("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) bt = 0.0;
   else
     bt = exp(gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
     return bt*betacf(a,b,x) / a;
  else
     return 1.0 - bt*betacf(b,a,1.0-x)/b;
}


double fcdf(double F, double v1, double v2) {

  return(1.0- betai(v2/2.0, v1/ 2.0, v2 / (v2 + v1* F)));
}

double studt(double t, double v) {
  return (1.0 - betai(v/(v+t*t),v/2.0, 0.5));
}

double pvnormal(double x) {
    // p-value of x ~ N(0,1)
    return 2.0*gaussianTail(fabs(x));
}

double pvbin(double p, int n, int k) {
 // return (1.0 - betai(k, n-k+1, p));
     return (1.0-betai(n-k, k+1, 1.0-p));
}


double gaussianTail(const double z)
{
  //  Based upon algorithm 5666 for the error function, from:
  //  Hart, J.F. et al, 'Computer Approximations', Wiley 1968
  static double p[7] = {220.2068679123761,
			221.2135961699311,
			112.0792914978709,
			33.91286607838300,
			6.373962203531650,
			0.7003830644436881,
			0.03526249659989109};
  static double q[8] = {440.4137358247522,
			793.8265125199484,
			637.3336333788311,
			296.5642487796737,
			86.78073220294608,
			16.06417757920695,
			1.755667163182642,
			0.08838834764831844};

  if(z > 37.)
    return 0.0;
  if(z < -37)
    return 1.0;

  double cutOff = 7.7071;
  double root2Pi = 2.506628274631001;
  double zabs = fabs(z);
  double expntl = exp(-.5*zabs*zabs);
  double pdf = expntl/root2Pi;
  double tail;
  if(zabs < cutOff)
    tail = expntl*((((((p[6]*zabs + p[5])*zabs + p[4])*zabs + p[3])*zabs + p[2])*zabs + p[1])*zabs + p[0])/(((((((q[7]*zabs + q[6])*zabs + q[5])*zabs + q[4])*zabs + q[3])*zabs + q[2])*zabs + q[1])*zabs + q[0]);

  else
    tail = pdf/(zabs + 1./(zabs + 2./(zabs + 3./(zabs + 4./(zabs + 0.65)))));

  return (z > 0.0) ? tail : 1.0-tail;

} // gaussianTail


