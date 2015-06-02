/*
gcc -Wall -ansi -pedantic /home/jsk/delt/src/fisi.c

R CMD SHLIB -o /home/jsk/fisi /home/jsk/delt/src/fisi.c

dyn.load("/home/jsk/fisi")
jako<-.C("findsplitCC",as.double(xdendat),
                       as.integer(obsoso),
                       as.integer(maara),
                       as.integer(curbeg),
                       as.integer(curend),
                       as.double(insuppo),
                       as.double(instep),
                       as.integer(ingrec),
                       as.integer(n),
                       as.integer(d),
                       as.integer(inmethod),
                       as.double(bound),
                       val = integer(1),
                       vec = integer(1),
                       #leftrec = integer(2*d+1),
                       #rightrec = integer(2*d+1),
                       leftbeg = integer(1),
                       leftend = integer(1),
                       rightbeg = integer(1),
                       rightend = integer(1),
                       obspoint = integer(maara+1))

apuu1 = double(1),apuu2 = double(1),apuu3 = double(1),
apuu4 = double(1))

*/

/*
double gapuu1, gapuu2, gapuu3, gapuu4;
*/

#include <math.h>
#include <stdlib.h> 


#define pieninfisi -2^1000000

double denssrCC(double, int, int, int);
int treesortbudCC(double *, int, int *);
int findobsCC(double *, int *, double, int);
int biggestindCC(int, double *);
double denvalCC(int, double, int);

void findsplitCC(double *xdendat,
                 int *obsoso,
                 int *maarain,
                 int *curbegin,
                 int *curendin,
                 double *suppo,
                 double *step,
                 int *grec, 
                 /*double *rec,*/
                 int *nin,
                 int *din,
                 int *methodin,
                 double *bound,
                 /* output */
		 int *valin,
                 int *vecin,
                 /*
                 int *gleftrec,
                 int *grightrec,
                 */
                 int *leftbegin,
                 int *leftendin,
                 int *rightbegin,
                 int *rightendin,
                 int *obspoint)

     /*
double *apuu1,
double *apuu2,
double *apuu3,
double *apuu4)
     */

{
    int maara, curbeg, curend, n, d, method;
    int i, j, k, l; 

    /*
    double xvec[*maarain+1];
    double ssrvec1[*din+1], ssrvec11[*din+1];  
                   ssrvec1:een talletetaan kullekin muuttujalle
                   pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                   jakopisteisiin liittyvista ssr arvoista) 
    double ssrvec2[*nin+2], ssrvec22[*nin+2];  
           matrix(1,jpislkm,1) 
	   ssrvec2:een talletet. kuhunkin mahd. 
	   jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
    double valvec[*din+1];    
    */
    double *xvec = (double *)malloc(sizeof(double) * (*maarain+1));
    double *ssrvec1 = (double *)malloc(sizeof(double) * (*din+1));
    double *ssrvec11 = (double *)malloc(sizeof(double) * (*din+1));

    double *ssrvec2 = (double *)malloc(sizeof(double) * (*nin+2));
    double *ssrvec22 = (double *)malloc(sizeof(double) * (*nin+2));
    double *valvec = (double *)malloc(sizeof(double) * (*din+1));

    /*
    int eligible1[*din+1];
    int dleftend[*din+1];             kullekin muuttujalle end of left 
    int gvalvec[*din+1];              kullekin muuttujalle jakopiste 
    int eligible2[*nin+2], elinum;
    int lefends[*nin+2];  matrix(1,jpislkm,1)
    int gleftrec[2*(*din)+1], grightrec[2*(*din)+1];
    int ordobsint[*maarain+1];  
    */
    int *eligible1 = (int *)malloc(sizeof(int) * (*din+1));
    int *dleftend = (int *)malloc(sizeof(int) * (*din+1));
    int *gvalvec = (int *)malloc(sizeof(int) * (*din+1));
    int *eligible2 = (int *)malloc(sizeof(int) * (*nin+2));
    int *lefends = (int *)malloc(sizeof(int) * (*nin+2));
    int *gleftrec = (int *)malloc(sizeof(int) * (2*(*din)+1));
    int *grightrec = (int *)malloc(sizeof(int) * (2*(*din)+1));
    int *ordobsint = (int *)malloc(sizeof(int) * (*maarain+1));

    int elinum;
    double suppvolume, supplen, valipit;
    double volumeleft, volumeright;
    int jpislkm;
    
    int rightbeg, rightend, leftbeg, leftend; 
    int gjakopiste;
    int minvali, apula, apu;
    int vec, gval; 
    int leftobslkm, rightobslkm;
    /* findobsCC */ 

    /*
    double fsdendat[*maarain+1][*din+1];
    int dordobs[*din+1][*maarain+1];   kullekin muuttujalle ordered obs 
    */
    double ** fsdendat;
    int ** dordobs;

    fsdendat = (double **)malloc((*maarain+1) * sizeof(double *));
    if (NULL == fsdendat) return;
    for (i = 0; i <= *maarain; i++) {
       fsdendat[i] = (double *)malloc((*din+1) * sizeof(double));
       if (NULL == fsdendat[i]) return;
    }

    dordobs = (int **)malloc((*din+1) * sizeof(int *));
    if (NULL == dordobs) return;
    for (i = 0; i <= *din; i++) {
       dordobs[i] = (int *)malloc((*maarain+1) * sizeof(int));
       if (NULL == dordobs[i]) return;
    }

    if ( xvec == NULL) return; 
    if ( ssrvec1 == NULL) return; 
    if ( ssrvec11 == NULL) return; 
    if ( valvec == NULL) return; 
    if ( eligible1 == NULL) return; 
    if ( dleftend == NULL) return; 
    if ( gvalvec == NULL) return; 
    if ( eligible2 == NULL) return; 
    if ( lefends == NULL) return; 
    if ( gleftrec == NULL) return; 
    if ( grightrec == NULL) return; 
    if ( ordobsint == NULL) return; 
    if ( ssrvec2 == NULL) return; 
    if ( ssrvec22 == NULL) return; 

    /* INITIALIZE */

    for (j=1; j<=*din; j++) eligible1[j]=1;
    for (j=1; j<=(*nin+1); j++) eligible2[j]=1;

    maara=*maarain;
    curbeg=*curbegin;
    curend=*curendin;
    n=*nin;
    d=*din;
    method=*methodin;
    for (j=1; j<=d; j++) gvalvec[j]=0;

    for (j=1; j<=maara; j++){
        for (l=1; l<=d; l++){
            fsdendat[j][l]=xdendat[(j-1)*d+l];
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

       jpislkm=grec[2*i]-grec[2*i-1]-1;    /*floor((n+1)*(valipit/supplen))-1;*/

       if (jpislkm>=1){   /* jos voidaan jakaa */

          for (j=1; j<=maara; j++){
              xvec[j]=fsdendat[j][i];
              ordobsint[j]=j;
          }

          apu=treesortbudCC(xvec,maara,ordobsint); 
          /* order pointers acc. to the i:th coord.*/
          /* treesortbudCC changes ordobsint (which is note book) */

          for (j=1; j<=maara; j++){
             dordobs[i][j]=ordobsint[j];
          }

          for (k=1; k<=jpislkm; k++){ 
                    /* kayd i:nnelle muuttujalle mahd. jakopist. lapi */

       	            /*jakopiste=supplen*k/(n+1);*/ 
	            gjakopiste=k;
 
        	    /* hila jaetaan vasempaan ja oikeaan hilaan */

	            /* leftrec=rec */
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

        	    /* kullekin jakopisteelle pointer to ordobs: */
                    /* last observation to belong to leftrec */

	            leftend=findobsCC(xvec,ordobsint,
                                      suppo[2*i-1]+gleftrec[2*i]*step[i],maara);
                    leftobslkm=leftend;
                    lefends[k]=leftend;
	            rightobslkm=maara-leftobslkm;

          	    /* volumeleft<-massone(leftrec) */
                    volumeleft=1;
                    for (j=1; j<=d; j++){
                       volumeleft=volumeleft*(gleftrec[2*j]-gleftrec[2*j-1])
                                  *step[j];
                    }

                    /* vas hilan estim arvo */
                    /*meanleft=denvalCC(leftobslkm,volumeleft,n);*/

         	    /* volumeright<-massone(rightrec) */
                    volumeright=1;
                    for (j=1; j<=d; j++){
                        volumeright=
                        volumeright*(grightrec[2*j]-grightrec[2*j-1])*step[j];
                    }

                    /* oik hilan estim arvo */
                    /*meanright=denvalCC(rightobslkm,volumeright,n);*/

                    ssrvec2[k]=denssrCC(volumeleft,leftobslkm,n,method)+
 	                       denssrCC(volumeright,rightobslkm,n,method);

		    if ((volumeleft/suppvolume<*bound) && (leftobslkm>0))
		      eligible2[k]=0;
		    if ((volumeright/suppvolume<*bound) && (rightobslkm>0))
		      eligible2[k]=0;

          }

          elinum=0;
          for (j=1; j<=jpislkm; j++){
	    if (eligible2[j]==1){
	       elinum=elinum+1;
	       ssrvec22[elinum]=ssrvec2[j];
	    }
          }
          /* haetaan indeksi, jossa ssr:n suurin arvo i. muuttuj. */
          if (elinum>0) minvali=biggestindCC(elinum,ssrvec22);
          else{
               minvali=biggestindCC(jpislkm,ssrvec2); 
               eligible1[j]=0;
          }
      
          gvalvec[i]=grec[2*i-1]+minvali; 
          ssrvec1[i]=ssrvec2[minvali];       /* min(ssrvec2) */
          dleftend[i]=lefends[minvali];

       }  /*(jpislkm>=1)*/ 

       else ssrvec1[i]=pieninfisi;

       i=i+1;
    }
 

  elinum=0;
  for (j=1; j<=d; j++){
      if (eligible1[j]==1){
	    elinum=elinum+1;
	    ssrvec11[elinum]=ssrvec1[j];
      }
  }
  /*haetaan indeksi, jossa ssr:n suurin arvo muuttujien yli*/
  if (elinum>0) vec=biggestindCC(elinum,ssrvec11);
  else vec=biggestindCC(d,ssrvec1); 
  /* sen muuttujan numero joka halkaistu */

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

 /* obspoint[beg:end]<-dordobs[vec,] */
 for (j=1; j<=maara; j++){
     apula=dordobs[vec][j];
     obspoint[j]=obsoso[apula];
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

     /*
     *apuu1=gapuu1;
     *apuu2=gapuu2;
     *apuu3=gapuu3;
     *apuu4=gapuu4;
     */

    free(xvec); 
    free(ssrvec1); 
    free(ssrvec11); 
    free(valvec);
    free(eligible1);
    free(dleftend);
    free(gvalvec);
    free(eligible2);
    free(lefends);
    free(gleftrec);
    free(grightrec);
    free(ordobsint);
    free(ssrvec2);
    free(ssrvec22);

    for(i = 0; i <= *maarain; i++) free(fsdendat[i]);
    free(fsdendat);
    for(i = 0; i <= *din; i++) free(dordobs[i]);
    free(dordobs);


}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double denssrCC(double volu, int nelem, int n, int method)

{
    double vastaus=0.0;

    if (nelem==0){
        vastaus=0;
    }
    else if (method==1){   /* 1 = loglik */
        vastaus=nelem*log(nelem/(n*volu));  /* log(N,base=exp) */
    }
    else{ /* if (method==2){  /* 2 = projec */  
        vastaus=pow(nelem,2)/(n*volu);
    }

return vastaus;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int treesortbudCC(double *vec, int numb, int *ordobsint)

{
/* ordobs is reference vector to be sorted*/

    int nodenum, i, node;
    int pinin, runner, apu, j;
    double curre;
    /*
    double value[2*numb], curre;
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

int findobsCC(double *ref, int *ordpointer, double leftrecend, int maara)

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

int biggestindCC(int maara, double *a)

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

double denvalCC(int maara, double masso, int n)

{
    double meano;

    meano=maara/(n*masso);

    return meano;
}




