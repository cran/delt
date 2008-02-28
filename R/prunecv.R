prunecv<-function(tree,pralfa,xtest,suppo,method){
#
# return vectors: for each alfa the nuber of alfas
# of the new subtrees which are on the interval,
# and the sum of estimation errors
# 
ntest<-dim(xtest)[1]
d<-dim(xtest)[2]
#
dime<-length(pralfa)
frek<-matrix(0,dime,1)
resu<-matrix(0,dime,1)
#
treeseq<-pruseq(tree)
#
subnum<-length(treeseq$alfa)
i<-1
while (i<=subnum){
  #
  curalfa<-treeseq$alfa[i]
  alloin<-lokeroi(pralfa,curalfa)
  #
  leafnum<-treeseq$leafs[i]
  puu<-picktree(treeseq,leafnum)
  #
  pv<-partition(puu,suppo)
  recs<-pv$recs
  values<-pv$values
  #
  recnum<-length(pv$values)
  #                    
  crit<-0
  j<-1
  while (j<=recnum){   #go through partition and calculate "crit"
     if (recnum>=2) rec<-recs[j,] else rec<-recs
     vol<-massone(rec)
     #                    
     # calculate the frequencies (of the test data) 
     #
     ob<-1
     fre<-0
     while (ob<=ntest){      #go through test set    
        obs<-xtest[ob,] 
        ans<-TRUE
        for (coordi in 1:d){
           if ((obs[coordi]<rec[2*coordi-1]) ||
	   (obs[coordi]>rec[2*coordi])) ans<-FALSE
        }
        if (ans) fre<-fre+1
        ob<-ob+1
     }
     # end of "fre" calculation
     if (method=="loglik"){
          if (fre>0) crit<-crit-fre*log(fre/(ntest*vol))/ntest
     }
     else{
          crit<-crit-fre^2/(ntest^2*vol)
     }
     j<-j+1
  }
  resu[alloin]<-resu[alloin]+crit
  frek[alloin]<-frek[alloin]+1
  i<-i+1
}
#
return(list(frek=frek,resu=resu))
}






