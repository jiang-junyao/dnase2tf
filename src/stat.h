/* 
 * File:   stat.h
 * Author: sjbaek
 *
 * Created on May 7, 2009, 2:06 PM
 */

#ifndef _STAT_H
#define	_STAT_H



#define MAXIT 100
#define EPS 1.0e-10
#define FPMIN 1.0e-30

double betacf(double a, double b, double x);
void nrerror(char* msg);
double gammln(double xx);
double betai(double a, double b, double x);
double fcdf(double F, double v1, double v2);
double studt(double t, double v);
double pvbin(double p, int n, int k);
double pvnormal(double x);
double gaussianTail(const double z);

#endif	/* _STAT_H */

