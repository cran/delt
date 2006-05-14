/*
gcc -Wall -ansi -pedantic /home/jsk/delt/src/muadhi.c


R CMD SHLIB -o /home/jsk/muadhi /home/jsk/delt/src/muadhi.c

dyn.load("/home/jsk/denfor/muadhi")
ds<-.C("densplitC",as.double(indendat),
                   as.integer(leaf),
                   as.integer(minlkm),
                   as.double(insuppo),
                   as.integer(inmethod),
                   as.integer(splitscan),
                   as.integer(seedf),
                   as.integer(n),
                   as.integer(d),
                   as.double(suppvol),
                   as.double(minvolume),
                   as.double(instep),
                   val = integer(maxnodnumR+1),
                   vec = integer(maxnodnumR+1),
                   mean = double(maxnodnumR+1),
                   nelem = integer(maxnodnumR+1),
                   ssr = double(maxnodnumR+1),
                   volume =  double(maxnodnumR+1),
                   left = integer(maxnodnumR+1),
                   right = integer(maxnodnumR+1),
                   glow = integer(d*maxnodnumR+1),
                   gupp = integer(d*maxnodnumR+1),
                   nodenum = integer(1),
                   obspointout = integer(maxnodnumR+1),
                   obslow = integer(maxnodnumR+1),
                   obsupp = integer(maxnodnumR+1))

apuu1 = double(1),apuu2 = double(1),apuu3 = double(1),
apuu4 = double(1))

*/

#include <math.h>
#include <stdlib.h> 

#define maxnodnum 10000
#define pienin -2^1000000
#define nmax 1000

/* GLOBALS */
int obspoint[nmax];   /* final ordering */
/* int ordobs[nmax];    note book for findsplitC */

/* output of findsplitC */
int gjakoval;
/*double jakoval; */
int jakovec;
int jakoleftbeg;
int jakoleftend;
int jakorightbeg;
int jakorightend;
int jakominvali;


double gapuu1;
double gapuu2;
double gapuu3;
double gapuu4;


/* declarations */
int findsplitC(double *,
              int *,
              int,
              int,
              int,
              double *,
              double *,
              int *,
	      /*double *,*/
              int ,
              int ,
              int);
double denssrC(double, int, int, int);
int treesortbudC(double *, int, int *);
int findobsC(double *, int *, double, int);
int biggestindC(int, double *);
double denvalC(int, double, int);


void densplitC(double *indendat,
               int *leaf,
               int *minlkm,
               double *suppo,
               int *method,
               int *splitscan,
               int *seedf,
               /* redundant */
               int *n,
               int *d, 
               double *suppvol,    /* volume of the support */  
               double *minvolume,  /* lowest resolution */ 
               double *step,
               /* output */
               int *val, 
               int *vec,
               double *mean,
               int *nelem,
               double *ssr,
               double *volume,
               int *left,
               int *right,
               int *glow,
               int *gupp,
               int *nodenum,
               int *obspointout,
               int *obslow,
               int *obsupp)

/*
double *apuu1,
double *apuu2,
double *apuu3,
double *apuu4)
*/
{

    /* int obspoint[*n+1];        pointers to the data matrix */
    int i, j, l, jin, pinin, curleafnum, curin, curparent;
    /*double currec[2*(*d)+1]; */
    int curbeg, curend, maara, obspointer;
    double masso;
    int jpistesti;
    double supplen, valipit;
    int jpislkm, test, jako;
    int rightbeg, rightend, leftbeg, leftend;
    /*double leftrec[2*(*d)+1], rightrec[2*(*d)+1];*/

    int pinoparen[maxnodnum+1];  /* stack of pointers to parents */
    int pinopointbeg[maxnodnum+1];
    int pinopointend[maxnodnum+1];
      /*pointers to obspoint, for each rectangle beg and and of the 
	locations of teh pointers to data*/

    /* 
    int obsoso[(*n)+1];
    int gcurrec[2*(*d)+1];
    int gleftrec[2*(*d)+1], grightrec[2*(*d)+1];
    */
    int *obsoso = (int *)malloc(sizeof(int) * (*n+1));
    int *gcurrec = (int *)malloc(sizeof(int) * (2*(*d)+1));
    int *gleftrec = (int *)malloc(sizeof(int) * (2*(*d)+1));
    int *grightrec = (int *)malloc(sizeof(int) * (2*(*d)+1));
    /*
    double xdendat[(*n)*(*d)+1];
    */
    double *xdendat = (double *)malloc(sizeof(double) * ((*n)*(*d)+1));
 
    /*
    double dendat[*n+1][*d+1];
    int ginternlow[maxnodnum+1][(*d)+1], ginternupp[maxnodnum+1][(*d)+1];
    int gpinorecs[maxnodnum+1][(*d)*2+1];  osioiden maaritelmat 
    */
    double ** dendat;
    int ** ginternlow;
    int ** ginternupp;
    int ** gpinorecs;

    dendat = (double **)malloc((*n+1) * sizeof(double *));
    if (NULL == dendat) exit(1);
    for (i = 0; i <= *n; i++) {
       dendat[i] = (double *)malloc((*d+1) * sizeof(double));
       if (NULL == dendat[i]) exit(1);
    }
    ginternlow = (int **)malloc((maxnodnum+1) * sizeof(int *));
    if (NULL == ginternlow) exit(1);
    for (i = 0; i <= maxnodnum; i++) {
       ginternlow[i] = (int *)malloc((*d+1) * sizeof(int));
       if (NULL == ginternlow[i]) exit(1);
    }
    ginternupp = (int **)malloc((maxnodnum+1) * sizeof(int *));
    if (NULL == ginternupp) exit(1);
    for (i = 0; i <= maxnodnum; i++) {
       ginternupp[i] = (int *)malloc((*d+1) * sizeof(int));
       if (NULL == ginternupp[i]) exit(1);
    }
    gpinorecs = (int **)malloc((maxnodnum+1) * sizeof(int *));
    if (NULL == gpinorecs) exit(1);    
    for (i = 0; i <= maxnodnum; i++) {
       gpinorecs[i] = (int *)malloc((2*(*d)+1) * sizeof(int));
       if (NULL == gpinorecs[i]) exit(1);
    }

    if ( obsoso == NULL) exit(1); 
    if ( gcurrec == NULL) exit(1); 
    if ( gleftrec == NULL) exit(1); 
    if ( grightrec == NULL) exit(1); 
    if ( xdendat == NULL) exit(1); 


    /* INITIALIZING */

    for (j=1; j<=(*n); j++){
        for (l=1; l<=(*d); l++){
            dendat[j][l]=indendat[(j-1)*(*d)+l];
        }
    }

    /* obspoint=seq(1:n) */
    /* all the observations belong to the root node */ 
    for (j=1; j<=*n; j++){
	obspoint[j]=j;
    }
    pinin=1;   /* pointer to the stacks */
    pinoparen[pinin]=0;
    /*
    for (j=1; j<=2*(*d); j++){
        pinorecs[pinin][j]=suppo[j];
    }
    */
    for (j=1; j<=(*d); j++){
        gpinorecs[pinin][2*j-1]=0;
        gpinorecs[pinin][2*j]=*n+1;
    }
    pinopointbeg[pinin]=1;
    pinopointend[pinin]=*n;

    curleafnum=1; /* current number of leafs */
    curin=0;      /* points to the tree-vectors to be created */ 

    /* MAIN LOOP */

while (pinin>=1){

    /* pull the data from the top of the stack to the result */

    curin=curin+1;
    curparent=pinoparen[pinin];
    for (j=1; j<=2*(*d); j++){
        gcurrec[j]=gpinorecs[pinin][j];
    }
    curbeg=pinopointbeg[pinin];
    curend=pinopointend[pinin];
    pinin=pinin-1;

    if (curparent>0){
        right[curparent]=curin;
    }
    val[curin]=0;  /* we have not splitted at any point */
    vec[curin]=0;  /* we have not splitted with respect to any direction */
    if (curbeg==0){ 
	maara=0;
    }
    else{
        maara=curend-curbeg+1; 
    }   
    nelem[curin]=maara;  /* the number of observations */

    masso=1;
    for (j=1; j<=*d; j++){
        masso=masso*(gcurrec[2*j]-gcurrec[2*j-1])*step[j];
    }

    volume[curin]=masso;   /* massone(currec) */

    mean[curin]=maara/((*n)*masso);
    /* value of the estimate for the rectangle */

    /*
    if (nelem[curin]==0){
        ssr[curin]=0;
    }
    else if (*method==1){ 
        ssr[curin]=nelem[curin]*log(nelem[curin]/(*n*volume[curin])); 
    }
    else if (*method==2){ 
        ssr[curin]=pow(nelem[curin],2)/(*n*volume[curin]);
    }
    */
    ssr[curin]=denssrC(volume[curin],nelem[curin],*n,*method);
    /* loglikeli of the estimate */

    obslow[curin]=curbeg;
    obsupp[curin]=curend;

    for (j=1; j<=*d; j++){
	ginternlow[curin][j]=gcurrec[2*j-1];
	ginternupp[curin][j]=gcurrec[2*j];
    }

    /* we go to the left subtree */

    jpistesti=0;   /* FALSE */
    for (j=1; j<=*d; j++){
	supplen=suppo[2*j]-suppo[2*j-1];
	   /* the length of the original interval */
        valipit=(gcurrec[2*j]-gcurrec[2*j-1])*step[j];
        jpislkm=floor((*n+1)*(valipit/supplen))-1;
        if (jpislkm>=1){
	    jpistesti=1;  /* TRUE */
        }
    }

    if ((*leaf)==0){
	test=((nelem[curin]>*minlkm) && (jpistesti));
    }
    else{
	test=((curleafnum<*leaf) && (nelem[curin]>*minlkm) && (jpistesti)); 
    }
    while (test==1){
 
        /* the node is to be splitted */
         
       curleafnum=curleafnum+1;
 
       for (j=curbeg; j<=curend; j++){
	   obspointer=obspoint[j];
           obsoso[j-curbeg+1]=obspointer;
           jin=j-curbeg+1;
           for (l=1; l<=(*d); l++){
               xdendat[(jin-1)*(*d)+l]=dendat[obspointer][l];
           }
       }

       if (*splitscan==0){
	   jako=findsplitC(xdendat,obsoso,
			   maara,curbeg,curend,suppo,step,gcurrec, /*currec,*/
                          *n,*d,*method);
       }
       else{
           jako=findsplitC(xdendat,obsoso,
			   maara,curbeg,curend,suppo,step,gcurrec, /*currec,*/
                          *n,*d,*method);
   
       }

       left[curin]=curin+1;

       val[curin]=gjakoval;
       vec[curin]=jakovec;
       leftbeg=jakoleftbeg;
       leftend=jakoleftend;
       rightbeg=jakorightbeg;
       rightend=jakorightend;

       /* leftrec<-rec */
       for (j=1; j<=2*(*d); j++){
           gleftrec[j]=gcurrec[j];
       }
       gleftrec[2*jakovec]=gjakoval;      /* vas rec:n loppupiste */

       /* rightrec<-rec */
       for (j=1; j<=2*(*d); j++){
           grightrec[j]=gcurrec[j];
       }
       grightrec[2*jakovec-1]=gjakoval;  /* oik rec:n alkupiste */

       /* the right child is pushed to the stack */

       pinin=pinin+1;
       pinoparen[pinin]=curin;
       for (j=1; j<=2*(*d); j++){
           gpinorecs[pinin][j]=grightrec[j];
       }
       pinopointbeg[pinin]=rightbeg;
       pinopointend[pinin]=rightend;
    
       /* left child is updated to the result */   

       curin=curin+1;

       val[curin]=0;  /* NA */
       vec[curin]=0;  /* NA */

       if (leftbeg==0){    /* 0 is a recognizer */ 
	  maara=0;
       }
       else{
          maara=leftend-leftbeg+1; 
       }   
       nelem[curin]=maara;  /* the number of observations */

       /* volume[curin]=massone(currec); */
       masso=1;
       for (j=1; j<=*d; j++){
            masso=masso*(gleftrec[2*j]-gleftrec[2*j-1])*step[j];
       }
       volume[curin]=masso;

       mean[curin]=denvalC(nelem[curin],volume[curin],*n);
       ssr[curin]=denssrC(volume[curin],nelem[curin],*n,*method);
     
       obslow[curin]=leftbeg;
       obsupp[curin]=leftend;

       for (j=1; j<=*d; j++){
           ginternlow[curin][j]=gleftrec[2*j-1];
           ginternupp[curin][j]=gleftrec[2*j];
       }

       /* for the possible further split of the left node */
       curbeg=leftbeg;
       curend=leftend;
       for (j=1; j<=2*(*d); j++){
           gcurrec[j]=gleftrec[j];
       }


       if (*leaf==0){
	   test=((nelem[curin]>*minlkm) && (volume[curin]>=*minvolume));
       } 
       else{
           test=((curleafnum<*leaf) &&
                 (nelem[curin]>*minlkm) && (volume[curin]>=*minvolume)); 
       }

    }  /* while (test) */

}  /* while (pinin>=1) */

 for (j=1; j<=curin; j++){
     for (l=1; l<=(*d); l++){
        glow[(j-1)*(*d)+l]=ginternlow[j][l];
        gupp[(j-1)*(*d)+l]=ginternupp[j][l];
     }
 }

 for (j=1; j<=*n; j++){
    obspointout[j]=obspoint[j];
 }

 *nodenum=curin;

/*
*apuu1=gapuu1;
*apuu2=gapuu2;
*apuu3=gapuu3;
*apuu4=gapuu4;
*/

    free(obsoso);
    free(gcurrec);
    free(gleftrec);
    free(grightrec);
    free(xdendat);

    for(i = 0; i <= *n; i++) free(dendat[i]);
    free(dendat);
    for(i = 0; i <= maxnodnum; i++) free(ginternlow[i]);
    free(ginternlow);
    for(i = 0; i <= maxnodnum; i++) free(ginternupp[i]);
    free(ginternupp);
    for(i = 0; i <= maxnodnum; i++) free(gpinorecs[i]);
    free(gpinorecs);

}


/*++++++++++++++++++++++++++++++++++++++++++++++++++*/


int findsplitC(double *xdendat,
              int *obsoso,
              int maara,
              int curbeg,
              int curend,
              double *suppo,
              double *step,
              int *grec,
	       /*double *rec,*/
              int n,
              int d,
              int method)

{
    int i, j, k, l;
    double suppvolume, supplen, valipit;
    double volumeleft, volumeright;
    int jpislkm;
    int rightbeg, rightend, leftbeg, leftend;
    double jakopiste;
    int gjakopiste;
    int minvali, vec, apu, apula, resu;
    double meanleft, meanright;
    int gval;
    int leftobslkm, rightobslkm;

    /*    
    double valvec[d+1];   
    double leftrec[2*d+1], rightrec[2*d+1]; 
    double xvec[maara+1];
    double ssrvec1[d+1];  ssrvec1:een talletetaan kullekin muuttujalle
                          pienin ssr-arvo (yli ko. muuttujan mahdollisiin
                          jakopisteisiin liittyvista ssr arvoista) 
    double ssrvec2[n+2];   matrix(1,jpislkm,1) 
	ssrvec2:een talletet. kuhunkin mahd. 
	jakopisteeseen liittyva ssr arvo i:nnelle muuttujalle 
    */
    double *valvec = (double *)malloc(sizeof(double) * (d+1));
    double *leftrec = (double *)malloc(sizeof(double) * (2*d+1));
    double *rightrec = (double *)malloc(sizeof(double) * (2*d+1));
    double *xvec = (double *)malloc(sizeof(double) * (maara+1));
    double *ssrvec1 = (double *)malloc(sizeof(double) * (d+1));
    double *ssrvec2 = (double *)malloc(sizeof(double) * (n+2));

    /*
    int dleftend[d+1];         kullekin muuttujalle end of left 
    int gvalvec[d+1];          kullekin muuttujalle jakopiste 
    int lefends[n+2];          matrix(1,jpislkm,1) 
    int ordobsint[maara+1];
    int gleftrec[2*d+1], grightrec[2*d+1];
    */
    int *dleftend = (int *)malloc(sizeof(int) * (d+1));
    int *gvalvec = (int *)malloc(sizeof(int) * (d+1));
    int *lefends = (int *)malloc(sizeof(int) * (n+2));
    int *ordobsint = (int *)malloc(sizeof(int) * (maara+1));
    int *gleftrec = (int *)malloc(sizeof(int) * (2*d+1));
    int *grightrec = (int *)malloc(sizeof(int) * (2*d+1));

    /*
    double fsdendat[maara+1][d+1];
    int dordobs[d+1][maara+1];  kullekin muuttujalle ordered obs 
    */
    double ** fsdendat;
    int ** dordobs;

    fsdendat = (double **)malloc((maara+1) * sizeof(double *));
    if (NULL == fsdendat) exit(1);
    for (i = 0; i <= maara; i++) {
       fsdendat[i] = (double *)malloc((d+1) * sizeof(double));
       if (NULL == fsdendat[i]) exit(1);
    }
    dordobs = (int **)malloc((d+1) * sizeof(int *));
    if (NULL == dordobs) exit(1);
    for (i = 0; i <= d; i++) {
       dordobs[i] = (int *)malloc((maara+1) * sizeof(int));
       if (NULL == dordobs[i]) exit(1);
    }

    if ( valvec == NULL) exit(1); 
    if ( leftrec == NULL) exit(1); 
    if ( rightrec == NULL) exit(1); 
    if ( xvec == NULL) exit(1); 
    if ( ssrvec1 == NULL) exit(1); 
    if ( ssrvec2 == NULL) exit(1);    
    if ( dleftend == NULL) exit(1); 
    if ( gvalvec == NULL) exit(1); 
    if ( lefends == NULL) exit(1); 
    if ( ordobsint == NULL) exit(1); 
    if ( gleftrec == NULL) exit(1); 
    if ( grightrec == NULL) exit(1); 
   


    /* INITIALIZE */

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
      supplen=suppo[2*i]-suppo[2*i-1];  /* alkup. jaettavan valin pituus */
      valipit=(grec[2*i]-grec[2*i-1])*step[i];
      /*valipit=rec[2*i]-rec[2*i-1];*/ 

      jpislkm=floor((n+1)*(valipit/supplen))-1;

      if (jpislkm>=1){   /* jos voidaan jakaa */

        for (j=1; j<=maara; j++){
            xvec[j]=fsdendat[j][i];
            ordobsint[j]=j;
        }

        apu=treesortbudC(xvec,maara,ordobsint); 
        /* order pointers acc. to the i:th coord.*/
        /* treesortbudC changes ordobsint (which is note book) */

        for (j=1; j<=maara; j++){
            dordobs[i][j]=ordobsint[j];
        }

        for (k=1; k<=jpislkm; k++){ 
            /* kayd i:nnelle muuttujalle mahd. jakopist. lapi */

	    jakopiste=supplen*k/(n+1); 
	    gjakopiste=k;

	    /* hila jaetaan vasempaan ja oikeaan hilaan */

	    /* leftrec=rec */
            for (j=1; j<=2*d; j++){
                /*leftrec[j]=rec[j];*/
                gleftrec[j]=grec[j];              
            }
            /*
	    leftrec[2*i-1]=rec[2*i-1];
	    leftrec[2*i]=rec[2*i-1]+jakopiste;
            */
	    gleftrec[2*i-1]=grec[2*i-1];
	    gleftrec[2*i]=grec[2*i-1]+gjakopiste; 

	    /* kullekin jakopisteelle pointer to ordobs: */
            /* last observation to belong to leftrec */

	    /*leftend=findobsC(xvec,ordobsint,leftrec[2*i],maara);*/
	    leftend=findobsC(xvec,ordobsint,
                             suppo[2*i-1]+gleftrec[2*i]*step[i],maara);
	    
            leftobslkm=leftend;
            lefends[k]=leftend;
	    rightobslkm=maara-leftobslkm;

            for (j=1; j<=2*d; j++){
                /*rightrec[j]=rec[j];*/
                grightrec[j]=grec[j];
            }
            /*
	    rightrec[2*i-1]=rec[2*i-1]+jakopiste;                  
	    rightrec[2*i]=rec[2*i];  
            */
	    grightrec[2*i-1]=grec[2*i-1]+gjakopiste;                  
	    grightrec[2*i]=grec[2*i];         
	    /* volumeleft<-massone(leftrec) */
            volumeleft=1;
            for (j=1; j<=d; j++){
                volumeleft=volumeleft*(gleftrec[2*j]-gleftrec[2*j-1])*step[j];
                /*volumeleft=volumeleft*(leftrec[2*j]-leftrec[2*j-1]);*/
            }

            /* vas hilan estim arvo */
            meanleft=denvalC(leftobslkm,volumeleft,n);

	    /* volumeright<-massone(rightrec) */
            volumeright=1;
            for (j=1; j<=d; j++){
             volumeright=volumeright*(grightrec[2*j]-grightrec[2*j-1])*step[j];
             /*volumeright=volumeright*(rightrec[2*j]-rightrec[2*j-1]);*/
            }

            /* oik hilan estim arvo */
            meanright=denvalC(rightobslkm,volumeright,n);

            ssrvec2[k]=denssrC(volumeleft,leftobslkm,n,method)+
 	               denssrC(volumeright,rightobslkm,n,method);

        }

	minvali=biggestindC(jpislkm,ssrvec2); 
        /* indeksi, jossa ssr:n suurin arvo i. muuttuj. */

        gvalvec[i]=grec[2*i-1]+minvali;
        /*valvec[i]=rec[2*i-1]+supplen*minvali/(n+1);*/
        ssrvec1[i]=ssrvec2[minvali];       /* min(ssrvec2) */
        dleftend[i]=lefends[minvali];

  }

  else ssrvec1[i]=pienin;

  i=i+1;
}
   
 vec=biggestindC(d,ssrvec1);
 /* sen muuttujan numero joka halkaistu */

 /*val=valvec[vec];*/
 gval=gvalvec[vec];        /* halkaisupiste */

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
     obspoint[curbeg+j-1]=obsoso[apula];
 }
     
/* 
resu<-list(val=val,vec=vec,leftrec=leftrec,rightrec=rightrec,leftbeg=leftbeg,
leftend=leftend,rightbeg=rightbeg,rightend=rightend,obspoint=obspoint)
*/

    gjakoval=gval;
    /*jakoval=val;*/ 
    jakovec=vec;
    jakoleftbeg=leftbeg;
    jakoleftend=leftend;
    jakorightbeg=rightbeg;
    jakorightend=rightend;
    jakominvali=minvali;

    resu=1;
    return resu;

    free( valvec );
    free( leftrec );
    free( rightrec );
    free( xvec );
    free( ssrvec1 );
    free( ssrvec2 );
    free( dleftend );
    free( gvalvec );
    free( lefends );
    free( ordobsint );
    free( gleftrec );
    free( grightrec );

    for(i = 0; i <= maara; i++) free(fsdendat[i]);
    free(fsdendat);
    for(i = 0; i <= d; i++) free(dordobs[i]);
    free(dordobs);

}

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

double denssrC(double volu, int nelem, int n, int method)

{
    double vastaus;

    if (nelem==0){
        vastaus=0;
    }
    else if (method==1){   /* 1 = loglik */
        vastaus=nelem*log(nelem/(n*volu));  /* log(N,base=exp) */
    }
    else if (method==2){  /* 2 = projec */
        vastaus=pow(nelem,2)/(n*volu);
    }

return vastaus;
}

/*++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int treesortbudC(double *vec, int numb, int *ordobsint)

{
/* ordobs is reference vector to be sorted*/

    int nodenum, i, node;
    double curre;
    int pinin, runner, apu, j;
    /*
    double value[2*numb], curre;
    int  left[2*numb], right[2*numb], refe[2*numb], pino[2*numb];    
    */
    double *value = (double *)malloc(sizeof(double) * (2*numb+1));
    int *left = (int *)malloc(sizeof(int) * (2*numb+1));
    int *right = (int *)malloc(sizeof(int) * (2*numb+1));
    int *refe = (int *)malloc(sizeof(int) * (2*numb+1));
    int *pino = (int *)malloc(sizeof(int) * (2*numb+1));

    if (value == NULL) exit(1); 
    if (left == NULL) exit(1); 
    if (right == NULL) exit(1); 
    if (refe == NULL) exit(1); 
    if (pino == NULL) exit(1); 


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

int findobsC(double *ref, int *ordpointer, double leftrecend, int maara)

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

  leftend=lcount;  /* could be zero */

return leftend;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

int biggestindC(int maara, double *a)

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

double denvalC(int maara, double masso, int n)

{
    double meano;

    meano=maara/(n*masso);

    return meano;
}


/*+++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* 
EI TARVITA 

int biggest_quicksort(void *a, int maara)

{

   int left, right, low, high;
   void *pivot_item; 
   double biggest, apu;

   low=1;
   high=maara;

   *pivot_item = a[low];
   pivot = left = low;
   right = high;
   while ( left < right ) {
     while( a[left] <= pivot_item ) left++;
     while( a[right] > pivot_item ) right--;
     if ( left < right ){
         apu=a[left];
         a[left]=a[right];
         a[right]=apu;
     }
   }
   a[low] = a[right];
   a[right] = pivot_item;

   biggest=a[maara];
   return biggest;
}
*/


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*

int makeorder(double *refe, int maara)

{

   int low, high, pivot, left, right;
   double pivot_item, apu;
   int appu, pivot_item_obs;   

   low=1;
   high=maara;

   pivot_item = refe[low];
   pivot = low;
   pivot_item_obs = ordobs[low];

   left = low;
   right = high;

   while ( left < right ){
     while( refe[left] <= pivot_item ) left++;
     while( refe[right] > pivot_item ) right--;
     if ( left < right ){
         apu=refe[left];
         refe[left]=refe[right];
         refe[right]=apu;

         appu=ordobs[left];
         ordobs[left]=ordobs[right];
         ordobs[right]=appu;
     }
   } 
   refe[low] = refe[right];
   refe[right] = pivot_item;
 
   ordobs[low] = ordobs[right];
   ordobs[right] = pivot_item_obs;

return right;
}

*/




