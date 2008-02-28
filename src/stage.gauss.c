/* 
gcc -Wall -ansi -pedantic /home/jsk/delt/src/stage.gauss.c


R CMD SHLIB -o /home/jsk/stage /home/jsk/delt/src/stage.gauss.c

dyn.load("/home/jsk/stage")

kg<-.C("stageGauss",
               as.double(mu0),
               as.double(sig0),
               as.double(indendat),
               as.integer(M),
               as.double(inmugrid),
               as.double(insiggrid),
               as.integer(insigeka),
               as.integer(insampstart),
               as.integer(n),
               as.integer(dictCard),
               as.integer(dictCardSig),
               as.integer(inboost),
               muut = double(M+1),
               sigit = double(M+1),
               curmix = double(M+1),

apu = double
)
kg$muut[2:(M+1)]
kg$sigit[2:(M+1)]
kg$curmix[2:(M+1)]

*/

#include <math.h>
#include <stdlib.h> 

#define maxn 2001
#define maxM 1001

double gaussprod(double mu1, double mu2, double sig1, double sig2);
double evanor(double x, int d);

void stageGauss(double *mu0,
                double *sig0,
                double *dendat,
                int    *M,
                double *mugrid,   
                double *siggrid,  
                int    *sigeka,
                int    *sampstart,
                int    *n,           /* length(dendat), sample size*/
                int    *dictCard,    /* length(mugrid) */
                int    *dictCardSig, /* length(siggrid) */
                int    *boost,
                /* output */
                double *muut,   /* means of the final mixture, M vec */
                double *sigit,  /* std:s of the final mixture, M vec */
                double *curmix  /* weights of the final mixture, M vec */
/*double *lapu */ 
)

{
 int i, j, ii, jj, k; 
 double piit[maxM+1];   /* matrix(0,M-1,1) */
 double sqint, val, point, evapoint;
 int imin, jmin, im, jm;
 double curmin;
 double prodint, gammanpik, gammanpikboost;
 double apu;
 /* 
 double riskit[*dictCard+1][*dictCardSig+1];
 */
 double ** riskit;
    riskit = (double **)malloc((*dictCard+1) * sizeof(double *));
    if (NULL == riskit) exit(1);
    for (i = 0; i <= *dictCard; i++) {
       riskit[i] = (double *)malloc((*dictCardSig+1) * sizeof(double));
       if (NULL == riskit[i]) exit(1);
    }

 /*apu = evanor(1,1); *lapu=apu;*/

 for (i=1; i<=(*M-1); i++) piit[i]=2.0/((double)i+2.0);

 if (*sampstart==1){
     imin=1;
     curmin=fabs(mugrid[1]-*mu0);
     for (i=1; i<=*dictCard; i++){
         apu=abs(mugrid[i]-*mu0);
         if (apu<curmin) imin=i;
     }
     jmin=1;
     curmin=abs(siggrid[1]-*sig0);
     for (j=1; j<=*dictCardSig; j++){
         apu=fabs(siggrid[j]-*sig0);
         if (apu<curmin) jmin=j;
     }
     muut[1]=mugrid[imin];
     sigit[1]=siggrid[jmin];
 }
 else{
    /* haetaan 1. termi */
    for (i=1; i<=*dictCard; i++){
       for (ii=1; ii<=*dictCardSig; ii++){
          sqint=gaussprod(0,0,siggrid[ii],siggrid[ii]);
          val=0;
          for (j=1; j<=*n; j++){   
	      point=dendat[j];
	      evapoint=(point-mugrid[i])/siggrid[ii];
	      val=val+evanor(evapoint,1)/siggrid[ii];
          }
          riskit[i][ii]=-2*val+sqint;
       }
    }

    imin=1;
    jmin=1;
    curmin=riskit[1][1];
    for (im=1; im<=*dictCard; im++){
        for (jm=1; jm<=*dictCardSig; jm++){
            if (riskit[im][jm]<=curmin){
	        curmin=riskit[im][jm]; 
                imin=im;
                jmin=jm;
            }
        }
    }
    /* mind=which.min(t(riskit))
       sarat=*dictCardSig;
       imin=ceiling(mind/sarat);
       jmin=mind-(imin-1)*sarat; */

    muut[1]=mugrid[imin];
    if (*sigeka==1) sigit[1]=1; else sigit[1]=siggrid[jmin]; 
 }

 /* haetaan termit 2-M */
 curmix[1]=1;             /* alussa vain yksi simppeli funktio */
 k=1;
 while (k<=(*M-1)){
   for (i=1; i<=*dictCard; i++){
      for (ii=1; ii<=*dictCardSig; ii++){
          /* calculate the -2*average of evaluations */
          val=0;
          for (j=1; j<=*n; j++){   
              point=dendat[j];
              evapoint=(point-mugrid[i])/siggrid[ii];
              val=val+evanor(evapoint,1)/siggrid[ii];
          }
        /* calculate the inner product of the candidate with the k-1 estimate*/
          prodint=0;
          jj=1;
          while (jj<=k){
	         apu=gaussprod(muut[jj],mugrid[i],sigit[jj],siggrid[ii]);
                 prodint=prodint+curmix[jj]*apu;
		 jj=jj+1;
          }
          /* calculate the risk at the k:th step */
          apu=gaussprod(0,0,siggrid[ii],siggrid[ii]);
          gammanpik=-2*piit[k]*val/(double)(*n)+pow(piit[k],2)*apu;
          gammanpikboost=-2*val/(double)(*n);
          if (*boost==0){
             riskit[i][ii]=gammanpik+2*(1-piit[k])*piit[k]*prodint;
          }    
          else{
            riskit[i][ii]=gammanpikboost+2*prodint;
	  } 
     }
   }  
   imin=1;
   jmin=1;
   curmin=riskit[1][1];
   for (im=1; im<=*dictCard; im++){
       for (jm=1; jm<=*dictCardSig; jm++){
           if (riskit[im][jm]<=curmin){
	       curmin=riskit[im][jm]; 
               imin=im;
               jmin=jm;
           }
       }
   }
   /* mind=which.min(t(riskit)) 
      sarat=*dictCardSig;
      imin=ceiling(mind/sarat);
      jmin=mind-(imin-1)*sarat; */

   muut[k+1]=mugrid[imin];
   sigit[k+1]=siggrid[jmin];

   for (i=1; i<=k; i++) curmix[i]=(1-piit[k])*curmix[i];
   /* curmix[1:k]=(1-piit[k])*curmix[1:k] */
   curmix[k+1]=piit[k];
   k=k+1;
 }   /*  while (k <= (M-1)) */


 for(i = 0; i <= *dictCard; i++) free(riskit[i]);
    free(riskit);

}  



double gaussprod(double mu1, double mu2, double sig1, double sig2)
{
    double result;

    result = pow(2*M_PI,-0.5)*pow(pow(sig1,2)+pow(sig2,2),-0.5)*
             exp(-pow(mu1-mu2,2)/(2*(pow(sig1,2)+pow(sig2,2))));

    return (result);
}


double evanor(double x, int d)
{
    double eta, normvakio, tulos;

    /* for (i=1; i<=d; i++) eta[i]<-pow(x[i],2);*/ 
    /*sum(x^2)=x:n pituuden nelio */
    eta=pow(x,2);
    normvakio=pow(sqrt(2*M_PI),-d); 
    tulos=exp(-0.5*eta)*normvakio;

    return tulos;
}
                           



