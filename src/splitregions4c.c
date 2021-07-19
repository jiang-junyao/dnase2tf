#include <R.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "stat.h"

#ifndef FALSE
#define FALSE (0)
#define TRUE (!(FALSE))
#endif

/// R CMD SHLIB splitregions4c.cpp stat.cpp


const int MAXFOOTPRINTS_HOTSPOT = 20000;
const int minHotspotWidth = 35;
const double trimRatio = 0.005;

double *ccount;
unsigned int *cmap;
int rstart; 
int rend;
			
			
//double fpbuffer[MAXFOOTPRINTS_HOTSPOT][4];


double mapfp(unsigned int a,unsigned int b ) 
{
  if (b<a || b<rstart || a>rend) 
     return 0.0;
       
  if (b>=rend)
      b = rend;

  if (b <rstart)
      return 0.0;

  if  (a<=rstart) 
            return cmap[b-rstart];
  else 
            return cmap[b-rstart] - cmap[a-rstart-1];
 
}

double sumfp(unsigned int a,unsigned int b)
{
  if (b<a || b<rstart || a>rend) 
     return 0.0;
       
  if (b>=rend)
      b = rend;

  if (b <rstart)
      return 0.0;

  if  (a<=rstart) 
            return ccount[b-rstart];
  else 
            return ccount[b-rstart] - ccount[a-rstart-1];
}


double getccount(unsigned int a) {
   if (a<rstart || a>rend) 
	return 0.0;
   else
        return ccount[a-rstart];
}


//#define getccount(a) (ccount[a-1])
//#define mapfp(a,b) ((double) (cmap[b-1]-cmap[a-2]))
    // Sum of cutcounts between a and b

//#define sumfp(a,b)  (ccount[b-1]-ccount[a-2])

void loadingtest() {
   printf("loading OK.\n");
}

void splitregions4c(double *ccount_arg,
    		unsigned int *cmap_arg,
 		int* prstart_arg, int* prend_arg,
   double* pprob, int* pminw, int* pmaxw,double* ppv_threshold,
   double* fpbuffer_st, 
   double* fpbuffer_ed,
   double* fpbuffer_cc,
   double* fpbuffer_pv,
   int* nfootprint) {


ccount = ccount_arg;
cmap = cmap_arg;
rstart = *prstart_arg;
rend = *prend_arg;

double prob = *pprob;
int minw = *pminw;
int maxw = *pmaxw;
double pv_threshold = *ppv_threshold;


    // Sum of cutcounts between a and b
//   cout << "findFootprint: rstart=" << rstart <<" rend=" << rend <<"\n";

    if (rend - rstart + 1 < minHotspotWidth) {
        //   plhs[0] = mxCreateDoubleMatrix(0,0,mxREAL);
  //      printf("rend-rstart+1<30\n");
        *nfootprint = 0;
        return;
    }

    int al, ar;
    int rstart2, rend2;
    double mincumCut, maxcumCut;
    int *dm, *dp;
    char *weaklm, *lm, *lp;

    double *diff, *z3, *active;

    int* Tst;
    int* Ted;
    double* Tcc;
    double* Tpv;
    int* Tprev;
    int* Tnext;
    char* Tneedcalc;

    int wst, wed;
    double lmap, wmap;
    double p1, p2, p3;

    int numregions;

    int mid1, mid2, mid3;
    int wst1, wst2, wst3;
    int wed1, wed2, wed3;
    double wmap1, wmap2, wmap3;
    double prob1, prob2, prob3;
    int *bl, *br;
    double n1, n2, n3;
    int curr, next;
    double numdots;

    bl = (int*) malloc((maxw + 1) * sizeof(int));
    br = (int*) malloc((maxw + 1) * sizeof(int));

//   mexPrintf("---->  CC=[%d,%d] mapsum=%d ccsum=%d\n",  
//          rstart,rend,  (unsigned int) mapfp(rstart,rend), (unsigned int) sumfp(rstart,rend) );

    al = (minw - 1) / 2;
    ar = minw / 2;
    for (int ii = 1; ii <= maxw; ii++) {
        bl[ii] = (int) (ii / prob * 0.5) - 1;
        br[ii] = bl[ii] + 1;
    }
    numdots = sumfp(rstart, rend);

    mincumCut = getccount(rstart)+ (numdots * trimRatio);
    maxcumCut = getccount(rend)- (numdots * trimRatio);

    rstart2 = rstart;
    for (int iii = rstart; iii <= rend; iii++) {
        if (getccount(iii) > mincumCut) {
            rstart2 = iii;
            break;
        }
    }
    rend2 = rend;
    for (int iii = rend; iii >= rstart; iii--) {
        if (getccount(iii) < maxcumCut) {
            rend2 = iii;
            break;
        }
    }

    //double W[MAXFOOTPRINTS_HOTSPOT+1][5];
    /// 1) Find "Seeds"
    int dlen;
    int *d;

    dlen = rend - rstart + 1;
    d = (int*) malloc((dlen + 1) * sizeof(int));
    int kk;
    if (rstart > 1) {
        for (int ii = rstart, kk=1; ii <= rend; ii++,kk++)
            d[kk] = (int) sumfp(ii, ii); //ccount[ii]-ccount[ii-1];
    }
    else {
        d[1] = (int) ccount[1];
        for (int ii = rstart+1,kk=2; ii <= rend; kk++,ii++)
            d[kk] = (int) sumfp(ii, ii); //ccount[ii] - ccount[ii-1];
    };

    ///%rstart2=
    ///% now d has the tag density

    dm = (int*) malloc((dlen + 1) * sizeof(int));
    dp = (int*) malloc((dlen + 1) * sizeof(int));
    weaklm = (char*) malloc((dlen + 1) * sizeof(char));


    lm = (char*) malloc((dlen + 1) * sizeof(char));
    lp = (char*) malloc((dlen + 1) * sizeof(char));


    dm[1] = 0;
    dp[dlen] = 0;

    for (int ii = 2; ii <= dlen - 1; ii++) {
        dm[ii] = d[ii - 1];
        dp[ii] = d[ii + 1];
        weaklm[ii] = (d[ii] <= dm[ii]) && (d[ii] <= dp[ii]);
        lm[ii] = (d[ii] < dm[ii]) & (d[ii] <= dp[ii]) || (d[ii] <= dm[ii]) & (d[ii] < dp[ii]);
    }
    //T=[];
    int mid, wsize;
    int st, ed;
    double N1, N2, NN;
    double z1, z2, pv;

    Tst = (int*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(int));
    Ted = (int*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(int));
    Tcc = (double*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(double)); // cutcounts
    Tpv = (double*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(double)); // p-values
    Tprev = (int*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(int)); // p-values
    Tnext = (int*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(int)); // p-values
    Tneedcalc = (char*) malloc((MAXFOOTPRINTS_HOTSPOT + 1) * sizeof(char));

    int jj;
    int ipos;
    ipos = rstart2;
    kk = 1;
    while (ipos <= rend2) {
        if (weaklm[ipos - rstart + 1]) {
            st = ipos;
            ed = st;
            for (jj = st; jj <= rend; jj++) {
                if ((lm[jj - rstart + 1] == 0 && weaklm[jj - rstart + 1] == 0) || (jj - st + 1 > maxw)) {
                    ed = jj - 1;
                    break;
                }
            }
            ipos = jj;

            mid = (st + ed) / 2;
            wsize = ed - st + 1;

            if (wsize > maxw) {
                ipos++;
                break;
            }

            wst = mid - bl[wsize];
            wed = mid + br[wsize];
            lmap = mapfp(st, ed);
            wmap = mapfp(wst, wed);

            if (wmap < 10 || lmap<1)
                prob3 = prob;
            else
                prob3 = lmap / wmap;

//            try {
                   n1= sumfp(st,ed);
                   N1= sumfp(wst, wed);
                   pv= (n1-N1*prob3)/sqrt(N1*prob3 * (1.0-prob3));
                  //  printf("(%d-%d) pv=%10.7e  wmap=%d, lmap=%d,  N=%d  n=%d  pr=%10.7e\n", wst,wed,pv,(int) wmap,(int) lmap,(int) N1,(int) n1, prob3);
  //          } catch (int e) {
  //              printf("error in binpvalue3\n"); // do nothing
   //         }
            //printf("(%d-%d) pv=%f\n", (int) wst, (int) wed, pv);
            //printf("(%d-%d) %10.7f\n", wst,wed,pv);
            Tst[kk] = st;
            Ted[kk] = ed;
            Tcc[kk] = sumfp(st, ed);
            Tpv[kk] = pv;

            Tprev[kk] = (kk == 1) ? -1 : (kk - 1);
            Tnext[kk] = kk + 1;
            Tneedcalc[kk] = TRUE;
            kk++;
            ipos++;
        } else {
            ipos++;
        }
    }

    numregions = kk - 1;

    if (numregions > 0) {
        Tnext[kk - 1] = -1;
        Tneedcalc[kk - 1] = FALSE;
    }


    //L= Ted-Tst+1;
    //[dummy, sorder] = sort(Ted);

    // 2) Merge

    //int L;
    //double ystart;
    int activeindex;


    activeindex = 0;
  //  ystart = -1 + 0.2;

    diff = (double*) malloc((numregions+ 1) * sizeof(double));
    z3 = (double*) malloc((numregions+ 1) * sizeof(double)) ;
    active = (double*) malloc((numregions+ 1) * sizeof(double));


    for (int ii = 1; ii <= numregions; ii++)
        active[ii] = Tnext[ii] > 0 || Tprev[ii] > 0;

    while (1) {

        int wsize1, wsize2, wsize3;

        for (int ii = 1; ii <= numregions; ii++) {
            if (!active[ii])
                diff[ii] = -1;
            else {
                activeindex = ii;
                curr = ii;
                next = Tnext[curr];
                if (Tneedcalc[curr] && next > 0) {
                    n1 = Tcc[curr];
                    wsize1 = Ted[curr] - Tst[curr] + 1;
                    n2 = Tcc[next];
                    wsize2 = Ted[next] - Tst[next] + 1;
                    n3 = sumfp(Tst[curr], Ted[next]);
                    wsize3 = Ted[next] - Tst[curr] + 1;

                    if (wsize3 > maxw) {
                        Tneedcalc[curr] = FALSE;
                        diff[curr] = -1.0;
                        continue;
                    }


                    mid1 = (Ted[curr] + Tst[curr]) / 2;
                    mid2 = (Ted[next] + Tst[next]) / 2;
                    mid3 = (Ted[next] + Tst[curr]) / 2;

                    wst1 = mid1 - bl[wsize1];
                    wed1 = mid1 + br[wsize1];
                    wst2 = mid2 - bl[wsize2];
                    wed2 = mid2 + br[wsize2];
                    wst3 = mid3 - bl[wsize3];
                    wed3 = mid3 + br[wsize3];

                    wmap1 = mapfp(wst1, wed1);
                    wmap2 = mapfp(wst2, wed2);
                    wmap3 = mapfp(wst3, wed3);

                    if (wmap3 < 5)
                        prob3 = prob;
                    else
                        prob3 = mapfp(Tst[curr], Ted[next]) / wmap3;

                    if (wmap2 < 5)
                        prob2 = prob;
                    else
                        prob2 = mapfp(Tst[next], Ted[next]) / wmap2;


                    if (wmap1 < 5)
                        prob1 = prob;
                    else
                        prob1 = mapfp(Tst[curr], Ted[curr]) / wmap1;


                   // try {
                        N1 = sumfp(mid1 - bl[wsize1], mid1 + br[wsize1]);
                        N2 = sumfp(mid2 - bl[wsize2], mid2 + br[wsize2]);
                        z1 = (n1 - N1 * prob1) / sqrt(N1 * prob1 * (1.0 - prob1));
                        z2 = (n2 - N2 * prob2) / sqrt(N2 * prob2 * (1.0 - prob2));

                        NN = sumfp(mid3 - bl[wsize3], mid3 + br[wsize3]);
                        if (wsize3 > maxw || NN < 10) {
                            z3[ii] = 1.0;
                            diff[ii] = 0;
                        } else {

                            z3[ii] = (n3 - NN * prob3) / sqrt(NN * prob3 * (1 - prob3));
                            p1 = (double) wsize1 / (double)wsize3;
                            p2 = (double) wsize2 / (double) wsize3;
                          //  p3 = 1.0 - p1 - p2;

                            diff[ii] = (p1 * z1 + p2 * z2 - (p1 + p2) * z3[ii]);
                            //%%%%diff[i]= (log2(wsize1/wsize3*z1 + wsize2/wsize3*z2)
                            //- log2((wsize1+wsize2)/wsize3*z3[i])) * (z3[i] < 0.2) *  (wsize3 <= maxw);

                        }
                   // }// try
                    //catch (int exct) {
                   //     printf("debug here\n");
                   // } // catch

                
              //       printf("(%d)-(%d) : diff=%f\n", (int) curr, (int) next, diff[curr]);
                
                    Tneedcalc[curr] = FALSE;
                } // if (Tneedcalc[curr] && next > 0)

            } // if (!active[ii]) ... else
        } // for (int=


        double dummy;
        int maxi;

        if (activeindex > 0) {
            diff[activeindex] = -1;
            Tneedcalc[activeindex] = FALSE;
            dummy = -1;
            for (int ll = 1; ll <= numregions; ll++) {
                if (dummy < diff[ll]) {
                    dummy = diff[ll];
                    maxi = ll;
                }
            }

          
          //  printf("[%d]-[%d] has the max. diff=%f\n", (int) maxi, (int) Tnext[maxi], diff[maxi]);
        } else dummy = 0.0;

        if (dummy > 0.0) {
            curr = maxi;
            next = Tnext[curr];

            //if next<0
            //	disp('debug here.');
            //end
            Tcc[curr] = sumfp(Tst[curr], Ted[next]);
            Tpv[curr] = z3[maxi];
            Ted[curr] = Ted[next];

            if (Tprev[curr] > 0 && Tnext[next] > 0) {
                Tnext[curr] = Tnext[next];
                Tprev[Tnext[next]] = curr;
            } else if (Tnext[next] < 0)
                Tnext[curr] = -1;
            else if (Tprev[curr] < 0) {
                Tnext[curr] = Tnext[next];
                Tprev[Tnext[next]] = curr;
            } // else if

            Tnext[next] = -1;
            Tprev[next] = -1;
            diff[next] = -1;
            Tneedcalc[curr] = TRUE;
            active[next] = FALSE;

            if (Tprev[curr] > 0)
                Tneedcalc[Tprev[curr]] = TRUE;
        } else
                break;  
    } // while 1


    // 3) Split


    int leng;
    double zscore;
    //double **footprint;

    kk = 1;
    for (int ii = 1; ii <= numregions; ii++) {
        if (Tnext[ii] > 0 || Tprev[ii] > 0) {
            leng = Ted[ii] - Tst[ii] + 1;
             zscore = Tpv[ii];

             if ((leng >= minw) & (zscore < pv_threshold)) {
                fpbuffer_st[kk-1] = Tst[ii];
                fpbuffer_ed[kk-1] = Ted[ii];
                fpbuffer_cc[kk-1] = Tcc[ii];
                fpbuffer_pv[kk-1] = Tpv[ii];
                //printf("%d: %d %d %f %f\n", kk, Tst[ii], Ted[ii], Tcc[ii], Tpv[ii]);
                kk = kk + 1;
            }
        }
    }

    *nfootprint = kk-1;
    free(Tst);
    free(Ted);
    free(Tcc);
    free(Tpv);
    free(Tprev);
    free(Tnext);
    free(Tneedcalc);
    free(diff);
    free(z3);
    free(active);

    
    free(bl);
    free(br);

    free(dm);
    free(dp);
    free(weaklm);
    free(lm);
    free(lp);
    
    return;
}
