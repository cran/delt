/*
gcc -Wall -ansi -pedantic /home/jsk/delt/src/fisipena.c

R CMD SHLIB -o /home/jsk/fisipena /home/jsk/delt/src/fisipena.c

dyn.load("/home/jsk/fisipena")
jako<-.C("findsplitpenaC",as.double(xdendat),
                          as.integer(obsoso),
                          as.integer(maara),
                          as.integer(curbeg),
                          as.integer(curend),
                          as.integer(n),
                          # mcdata
                          as.double(mcxdendat),
                          as.integer(mcobsoso),
                          as.integer(mcmaara),
                          as.integer(mccurbeg),
                          as.integer(mccurend),
                          as.integer(mcn),
                          # general
                          as.double(insuppo),
                          as.double(instep),
                          as.integer(ingrec),
                          as.integer(d),
                          as.integer(inmethod),
                          as.double(bound),
                          as.double(mix),
                          # output
                          val = integer(1),
                          vec = integer(1),
                          #
                          leftbeg = integer(1),
                          leftend = integer(1),
                          rightbeg = integer(1),
                          rightend = integer(1),
                          obspoint = integer(maara+1),
                          # same for mc-sample 
                          mcleftbeg = integer(1),
                          mcleftend = integer(1),
                          mcrightbeg = integer(1),
                          mcrightend = integer(1),
                          mcobspoint = integer(maara+1))

apuu1 = double(1),apuu2 = double(1),apuu3 = double(1),
apuu4 = double(1))

*/

/*
double gapuu1, gapuu2, gapuu3, gapuu4;
*/

#include <math.h>
#include <stdlib.h> 

#define pieninpena -2^1000000

double denssrCCC(double, int, int, int, double);
int treesortbudCCC(double *, int, int *);
int findobsCCC(double *, int *, double, int);
int biggestindCCC(int, double *);
double denvalCCC(int, double, int);

void findsplitpenaC(double *xdendat,
                    int *obsoso,
                    int *maarain,                 
                    int *curbegin,
                    int *curendin,
                    int *nin,
                    /*mcdata*/
		    double *mcxdendat,
                    int *mcobsoso,
                    int *mcmaarain,                 
                    int *mccurbegin,
                    int *mccurendin, 
                    int *mcnin,
		    /*general*/
                    double *suppo,
                    double *step,
                    int *grec, 
                    int *din,
                    int *methodin,
                    double *bound,
                    double *mixin,

                    /* output */
		    int *valin,
                    int *vecin,
                    
                    int *leftbegin,
                    int *leftendin,
                    int *rightbegin,
                    int *rightendin,
                    int *obspoint,
                    /* same for mc-sample */
                    int *mcleftbegin,
                    int *mcleftendin,
                    int *mcrightbegin,
                    int *mcrightendin,
                    int *mcobspoint,
double *apuu1,
double *apuu2,
double *apuu3,
double *apuu4)


{

    int d, method;
    int i, j, k, l;

    /* these are needed both original and mc-sample */
    int maara, curbeg, curend, n;
    int mcmaara, mccurbeg, mccurend, mcn;
    /* kullekin muuttujalle ordered obs */
    int leftobslkm, rightobslkm, mcleftobslkm, mcrightobslkm;
    int rightbeg, rightend, leftbeg, leftend; 
    int mcrightbeg, mcrightend, mcleftbeg, mcleftend; 
   
    double suppvolume, supplen, valipit;
    double volumeleft, volumeright;
    int jpislkm;
    int gjakopiste;
    int minvali, apu, apula;
    int vec, gval; 
    double meanleft, meanright;
    /* findobsCCC */ 
    double mix;

    /*
    double xvec[*maarain+1];
    double  mcxvec[*mcmaarain+1];
       
    double ssrvec1[*din+1], ssrvec11[*din+1];  
    double ssrvec2[*nin+2], ssrvec22[*nin+2];  
    */
       /*
       ssrvec1:een talletetaan kullekin muuttujalle
       pienin ssr-arvo (yli ko. muuttujan mahdollisiin
       jakopisteisiin liittyvista ssr arvoista) 

       ssrvec2:een talletet. kuhunkin mahd. 
       jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
       */
    
    double *xvec = (double *)malloc(sizeof(double) * (*maarain+1));
    
    double *mcxvec = (double *)malloc(sizeof(double) * (*mcmaarain+1));
       
    double *ssrvec1 = (double *)malloc(sizeof(double) * (*din+1));
    double *ssrvec11 = (double *)malloc(sizeof(double) * (*din+1));
    double *ssrvec2 = (double *)malloc(sizeof(double) * (*nin+2));
    double *ssrvec22 = (double *)malloc(sizeof(double) * (*nin+2));

    /*
    int eligible1[*din+1];
    int eligible2[*nin+2];
   
    int dleftend[*din+1], mcdleftend[*din+1];
    int gvalvec[*din+1];           
   
    int lefendsit[*nin+2], mcmclefends[*nin+2];
    int gleftrec[2*(*din)+1], grightrec[2*(*din)+1]; 
    */ 
    /*
    int ordobsint[*maarain+1];
    int  mcordobsint[*mcmaarain+1]; 
    */
    int *eligible1 = (int *)malloc(sizeof(int) * (*din+1));
    int *eligible2 = (int *)malloc(sizeof(int) * (*nin+2));
    int *dleftend = (int *)malloc(sizeof(int) * (*din+1));
    int *mcdleftend = (int *)malloc(sizeof(int) * (*din+1));
    /*  kullekin muuttujalle end of left */  
    int *gvalvec = (int *)malloc(sizeof(int) * (*din+1));
    /* kullekin muuttujalle jakopiste */  
    int *lefendsit = (int *)malloc(sizeof(int) * (*nin+2));
    int *mcmclefends = (int *)malloc(sizeof(int) * (*nin+2));
    int *gleftrec = (int *)malloc(sizeof(int) * (2*(*din)+1));
    int *grightrec = (int *)malloc(sizeof(int) * (2*(*din)+1));
    
     
    int *ordobsint = (int *)malloc(sizeof(int) * (*maarain+1));
    
    int *mcordobsint = (int *)malloc(sizeof(int) * (*mcmaarain+1));
    

    /*
    double fsdendat[*maarain+1][*din+1], mcfsdendat[*mcmaarain+1][*din+1];
    int dordobs[*din+1][*maarain+1], mcdordobs[*din+1][*mcmaarain+1];  
    */
    double ** fsdendat;
    double ** mcfsdendat;
    int ** dordobs;
    int ** mcdordobs;

    fsdendat = (double **)malloc((*maarain+1) * sizeof(double *));
    if (NULL == fsdendat) return;
    for (i = 0; i <= *maarain; i++) {
       fsdendat[i] = (double *)malloc((*din+1) * sizeof(double));
       if (NULL == fsdendat[i]) return;
    }

    mcfsdendat = (double **)malloc((*mcmaarain+1) * sizeof(double *));
    if (NULL == mcfsdendat) return;
    for (i = 0; i <= *mcmaarain; i++) {
       mcfsdendat[i] = (double *)malloc((*din+1) * sizeof(double));
       if (NULL == mcfsdendat[i]) return;
    }

    dordobs = (int **)malloc((*din+1) * sizeof(int *));
    if (NULL == dordobs) return;
    for (i = 0; i <= *din; i++) {
       dordobs[i] = (int *)malloc((*maarain+1) * sizeof(int));
       if (NULL == dordobs[i]) return;
    }

    mcdordobs = (int **)malloc((*din+1) * sizeof(int *));
    if (NULL == mcdordobs) return;
    for (i = 0; i <= *din; i++) {
       mcdordobs[i] = (int *)malloc((*mcmaarain+1) * sizeof(int));
       if (NULL == mcdordobs[i]) return;
    }

    
    if ( xvec == NULL) return; 
    
    if ( mcxvec == NULL) return; 
    
    if ( ssrvec1 == NULL) return; 
    if ( ssrvec11 == NULL) return; 
    if ( ssrvec2 == NULL) return; 
    if ( ssrvec22 == NULL) return; 

    
    if ( eligible1 == NULL) return;
    if ( eligible2 == NULL) return; 
    
     
    if ( dleftend == NULL) return; 
    if ( mcdleftend == NULL) return; 
    if ( gvalvec == NULL) return; 
    
    
    if ( lefendsit == NULL) return; 
    if ( mcmclefends == NULL) return; 
    if ( gleftrec == NULL) return; 
    if ( grightrec == NULL) return; 
    
    
    if ( ordobsint == NULL) return; 
    
    if ( mcordobsint == NULL) return; 
    

    /* INITIALIZE */

    for (j=1; j<=*din; j++) eligible1[j]=1;
    for (j=1; j<=(*nin+1); j++) eligible2[j]=1;

    mix=*mixin;
    maara=*maarain;
    curbeg=*curbegin;
    curend=*curendin;
    n=*nin;

    mcmaara=*mcmaarain;
    mccurbeg=*mccurbegin;
    mccurend=*mccurendin;
    mcn=*mcnin; 

    d=*din;
    method=*methodin;
    for (j=1; j<=d; j++) gvalvec[j]=0;

    for (j=1; j<=maara; j++){
        for (l=1; l<=d; l++){
            fsdendat[j][l]=xdendat[(j-1)*d+l];
        }
    }
    for (j=1; j<=mcmaara; j++){
        for (l=1; l<=d; l++){
            mcfsdendat[j][l]=mcxdendat[(j-1)*d+l];
        }
    }

    /* suppvolume<-massone(suppo) */
    suppvolume=1;
    for (j=1; j<=d; j++){
        suppvolume=suppvolume*(suppo[2*j]-suppo[2*j-1]);
    }

    i=1;
    while (i<=d){     /* kaydaan muuttujat lapi */

       supplen=suppo[2*i]-suppo[2*i-1];/*alkup.jaettavan valin pit*/
       valipit=(grec[2*i]-grec[2*i-1])*step[i]; 

       jpislkm=grec[2*i]-grec[2*i-1]-1;  /*floor((n+1)*(valipit/supplen))-1;*/

       if (jpislkm>=1){   /* jos voidaan jakaa */

          for (j=1; j<=jpislkm; j++) lefendsit[j]=0;

          for (j=1; j<=maara; j++){
              xvec[j]=fsdendat[j][i];
              ordobsint[j]=j;
          }
          /* same for mc-sample */
          for (j=1; j<=mcmaara; j++){
              mcxvec[j]=mcfsdendat[j][i];
              mcordobsint[j]=j;
          }

          apu=treesortbudCCC(xvec,maara,ordobsint); 
          /* order pointers acc. to the i:th coord.*/
          /* treesortbudCCC changes ordobsint (which is note book) */
          for (j=1; j<=maara; j++){
             dordobs[i][j]=ordobsint[j];
          }
          /* same for mc-sample */
          apu=treesortbudCCC(mcxvec,mcmaara,mcordobsint); 
          for (j=1; j<=mcmaara; j++){
             mcdordobs[i][j]=mcordobsint[j];
          }
        
          k=1;
          while (k<=jpislkm){
          /*for (k=1; k<=jpislkm; k++){ */
                    /* kayd i:nnelle muuttujalle mahd. jakopist. lapi */

	            gjakopiste=k;
 
        	    /* hila jaetaan vasempaan ja oikeaan hilaan */

                    for (j=1; j<=2*d; j++){
                        gleftrec[j]=grec[j];              
                    }
	            gleftrec[2*i-1]=grec[2*i-1];
	            gleftrec[2*i]=grec[2*i-1]+gjakopiste; 

                    for (j=1; j<=2*d; j++){
                       grightrec[j]=grec[j]; 
                    }
	            grightrec[2*i-1]=grec[2*i-1]+gjakopiste;                  
	            grightrec[2*i]=grec[2*i];         

          	    /* volumeleft<-massone(leftrec) */
                    volumeleft=1;
                    for (j=1; j<=d; j++){
                      volumeleft=volumeleft*(gleftrec[2*j]-gleftrec[2*j-1])
                                 *step[j];
                    }
         	    /* volumeright<-massone(rightrec) */
                    volumeright=1;
                    for (j=1; j<=d; j++){
                        volumeright=
                            volumeright*(grightrec[2*j]-grightrec[2*j-1])*step[j];
                    }

        	    /* kullekin jakopisteelle pointer to ordobs: */
                    /* last observation to belong to leftrec */
	            leftend=findobsCCC(xvec,ordobsint,
                                       suppo[2*i-1]+gleftrec[2*i]*step[i],maara);
	            leftobslkm=leftend;
                    lefendsit[gjakopiste]=leftend;
	            rightobslkm=maara-leftobslkm;

                    /* same for mc-sample */
	            mcleftend=findobsCCC(mcxvec,mcordobsint,
                              suppo[2*i-1]+gleftrec[2*i]*step[i],mcmaara);
	            mcleftobslkm=mcleftend;
                    mcmclefends[k]=mcleftend;
	            mcrightobslkm=mcmaara-mcleftobslkm;

                    /* vas hilan estim arvo */
                    meanleft=denvalCCC(leftobslkm,volumeleft,n);
                    /* oik hilan estim arvo */
                    meanright=denvalCCC(rightobslkm,volumeright,n);


                    ssrvec2[k]=denssrCCC(volumeleft,leftobslkm,n,method,mix)+denssrCCC(volumeright,rightobslkm,n,method,mix)+2*(mix-1)*mix*mcleftobslkm*meanleft/mcn+2*(mix-1)*mix*mcrightobslkm*meanright/mcn;
 
                    /*
		    if ((volumeleft/suppvolume<*bound) && (leftobslkm>0))
		      eligible2[k]=0;
		    if ((volumeright/suppvolume<*bound) && (rightobslkm>0))
		      eligible2[k]=0;
		    */

		    k=k+1;
          }  /* while k<= jpislkm */

          /*
          elinum=0;
          for (j=1; j<=jpislkm; j++){
	    if (eligible2[j]==1){
	      elinum=elinum+1;
	      ssrvec22[elinum]=ssrvec2[j];
	    }
          }
          */

          /* haetaan indeksi, jossa ssr:n suurin arvo i. muuttuj. */
          /* 
         if (elinum>0) minvali=biggestindCCC(elinum,ssrvec22);
          else{
               minvali=biggestindCCC(jpislkm,ssrvec2); 
               eligible1[j]=0;
          }
          */
          minvali=biggestindCCC(jpislkm,ssrvec2); 

          gvalvec[i]=grec[2*i-1]+minvali; 
          ssrvec1[i]=ssrvec2[minvali];       /* min(ssrvec2) */

          dleftend[i]=lefendsit[minvali];
          mcdleftend[i]=mcmclefends[minvali];

       }  /*if (jpislkm>=1)*/ 

       else ssrvec1[i]=pieninpena;

       i=i+1;
    }
   
    /*
  elinum=0;
  for (j=1; j<=d; j++){
      if (eligible1[j]==1){
	    elinum=elinum+1;
	    ssrvec11[elinum]=ssrvec1[j];
      }
  }
    */

  /*haetaan indeksi, jossa ssr:n suurin arvo muuttujien yli*/
  /*
  if (elinum>0) vec=biggestindCCC(elinum,ssrvec11);
  else vec=biggestindCCC(d,ssrvec1); 
  */
  /* sen muuttujan numero joka halkaistu */
  vec=biggestindCCC(d,ssrvec1); 

 /*val=valvec[vec];*/ 
  gval=gvalvec[vec];        /*halkaisupiste */

  if (dleftend[vec]==0){ /* no observations in the left rec */
     leftbeg=0;
     leftend=0;
     rightbeg=curbeg;
     rightend=curend;
  }
  else if (dleftend[vec]==maara){ /* no observations in the right rec */
     leftbeg=curbeg;
     leftend=curend;
     rightbeg=0;
     rightend=0;
  }
  else{
      leftbeg=curbeg;
      leftend=curbeg+dleftend[vec]-1;
      rightbeg=leftend+1;
      rightend=curend;  
  }

  /* same for mc-sample */
  if (mcdleftend[vec]==0){ /* no observations in the left rec */
     mcleftbeg=0;
     mcleftend=0;
     mcrightbeg=mccurbeg;
     mcrightend=mccurend;
  }
  else if (mcdleftend[vec]==mcmaara){ /* no observations in the right rec */
     mcleftbeg=mccurbeg;
     mcleftend=mccurend;
     mcrightbeg=0;
     mcrightend=0;
  }
  else{
      mcleftbeg=mccurbeg;
      mcleftend=mccurbeg+mcdleftend[vec]-1;
      mcrightbeg=mcleftend+1;
      mcrightend=mccurend;  
  }


 /* obspoint[beg:end]<-dordobs[vec,] */
 for (j=1; j<=maara; j++){
     apula=dordobs[vec][j];
     obspoint[j]=obsoso[apula];
 }

 for (j=1; j<=mcmaara; j++){
     apula=mcdordobs[vec][j];
     mcobspoint[j]=mcobsoso[apula];
 }
   
/* 
resu<-list(val=val,vec=vec,leftrec=leftrec,rightrec=rightrec,leftbeg=leftbeg,
leftend=leftend,rightbeg=rightbeg,rightend=rightend,obspoint=obspoint)
*/

     *valin=gval;
     *vecin=vec;
     *leftbegin=leftbeg;
     *leftendin=leftend;
     *rightbegin=rightbeg;
     *rightendin=rightend;

     *mcleftbegin=mcleftbeg;
     *mcleftendin=mcleftend;
     *mcrightbegin=mcrightbeg;
     *mcrightendin=mcrightend;

/*     
     *apuu1=gapuu1;
     *apuu2=gapuu2;
     *apuu3=gapuu3;
     *apuu4=gapuu4;
     */

    for(i = 0; i <= *maarain; i++) free(fsdendat[i]);
    free(fsdendat);
    for(i = 0; i <= *mcmaarain; i++) free(mcfsdendat[i]);
    free(mcfsdendat);
    for(i = 0; i <= *din; i++) free(dordobs[i]);
    free(dordobs);
    for(i = 0; i <= *din; i++) free(mcdordobs[i]);
    free(mcdordobs);

    
    free(xvec); 
    
    free(mcxvec); 
     
    free(ssrvec1); 
    free(ssrvec11); 
    free(ssrvec2);
    free(ssrvec22);

        
    free(eligible1);
    free(eligible2);
    
    
    free(dleftend);
    free(mcdleftend);
    free(gvalvec);
    
    
    free(lefendsit);
    free(mcmclefends);
    free(gleftrec);
    free(grightrec);
    
    
    free(ordobsint);
    
    free(mcordobsint);
    
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double denssrCCC(double volu, int nelem, int n, int method, double mix)

{
    double vastaus=0.0;

    if (nelem==0){
        vastaus=0;
    }
    else if (method==1){   /* 1 = loglik */
        vastaus=nelem*log(nelem/(n*volu))/n;  /* log(N,base=exp) */
    }
    else{ /* if (method==2){  /* 2 = projec */
        vastaus=mix*(2-mix)*pow(nelem,2)/(pow(n,2)*volu);
    }

    return vastaus;

}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int treesortbudCCC(double *vec, int numb, int *ordobsint)

{
/* ordobs is reference vector to be sorted*/

    int nodenum, i, node;
    int pinin, runner, apu, j;
    double curre;
    /*
    double value[2*numb];
    int  left[2*numb], right[2*numb], refe[2*numb], pino[2*numb];    
    */
    double *value = (double *)malloc(sizeof(double) * (2*numb+1));
    int *left = (int *)malloc(sizeof(int) * (2*numb+1));
    int *right = (int *)malloc(sizeof(int) * (2*numb+1));
    int *refe = (int *)malloc(sizeof(int) * (2*numb+1));
    int *pino = (int *)malloc(sizeof(int) * (2*numb+1));

    if (value == NULL) return 0; 
    if (left == NULL) return 0; 
    if (right == NULL) return 0; 
    if (refe == NULL) return 0; 
    if (pino == NULL) return 0; 


    /* initialize */
    for (j=1; j<=(2*numb); j++){
	left[j]=0;
    }

    value[1]=vec[1];
    refe[1]=1;
    nodenum=1;

    i=2;
    while (i<=numb){
  
	curre=vec[i]; 
	node=1;

	/* go to a leaf */
        while (left[node]>0){
            if (curre<value[node]){
                   node=left[node];
            } 
            else{ 
		node=right[node];
            }
        }

	/* where are in a leaf: add 2 new leafs */
        if (curre<value[node]){  
	    value[nodenum+1]=curre;
	    value[nodenum+2]=value[node];

            refe[nodenum+1]=i;
            refe[nodenum+2]=refe[node];

        }
        else{
	    value[nodenum+1]=value[node];
	    value[nodenum+2]=curre;

            refe[nodenum+1]=refe[node];
            refe[nodenum+2]=i;

        }

	left[node]=nodenum+1;
	right[node]=nodenum+2;
	value[node]=(curre+value[node])/2;

	nodenum=nodenum+2;

	i=i+1;      
    }

    /* collect sorted vector */

    pino[1]=1;
    pinin=1;
    runner=1;

    while (pinin>0){

	/* take from stack */

	node=pino[pinin];
	pinin=pinin-1;

	if (left[node]==0){  /* we are in a leaf */
	    ordobsint[runner]=refe[node];
	    runner=runner+1;
        }

        while (left[node]>0){

	    pinin=pinin+1;
	    pino[pinin]=right[node];

	    node=left[node];

	    if (left[node]==0){  /* we are in a leaf */
		ordobsint[runner]=refe[node];
		runner=runner+1;
            }
 
        }
 
    }

    apu=1;
    return apu;


    free(value);
    free(left);
    free(right);
    free(refe);
    free(pino);


}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++*/

int findobsCCC(double *ref, int *ordpointer, double leftrecend, int maara)

{
  int lcount, leftend, j, ind;
  double havakoor;

  lcount=0;
  j=1;
  ind=ordpointer[j];
  havakoor=ref[ind];
  
  while ((havakoor<leftrecend) && (j<maara)){
        lcount=lcount+1;

        j=j+1;
        ind=ordpointer[j];
        havakoor=ref[ind];
  }

  ind=ordpointer[maara];
  havakoor=ref[ind];
  if (havakoor<=leftrecend) lcount=lcount+1;

  leftend=lcount;  /* could be zero */

return leftend;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int biggestindCCC(int maara, double *a)

{
    int j, bigind;
    double biggest, curre;

    biggest=a[1];
    bigind=1;
    for (j=1; j<=maara; j++){
        curre=a[j];
        if (curre>biggest){
	    biggest=curre;
            bigind=j;
        }
    }     

    return bigind;
}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double denvalCCC(int maara, double masso, int n)

{
    double meano;

    meano=maara/(n*masso);

    return meano;
}




