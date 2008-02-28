crosvali<-function(wv,dendat,minlkm,suppo,seed,method="loglik",blokki=50){
#cartcv.xpl
#wv: integer >=2, wv fold cross-validation is performed, that
#              is, the data is divided in wv number of ways to an 
#              estimation set and a test set.
#              Division is formed randomly.

n<-dim(dendat)[1]
ositteet<-osita(n,wv,seed)
tr<-densplit(dendat,minlkm,suppo,method=method)  #,blokki)

prseq<-pruseq(tr)
pralfa<-prseq$alfa
lnumber<-prseq$leafs
tuldim<-length(pralfa)

cv<-matrix(0,tuldim,1)
cvstd<-matrix(0,tuldim,1)
numinmean<-matrix(0,tuldim,1)

apu<-matrix(NA,wv,tuldim)
i<-1
while (i<=wv){
  osite<-ositalog(n,ositteet,i)
  kaanto<-matrix(1,n,1)-osite
  xdata<-paf(dendat,kaanto)
  suppo<-supp(xdata)
  tree<-densplit(xdata,minlkm,suppo,method=method)  #,blokki)
  xtest<-paf(dendat,osite)
  pc<-prunecv(tree,pralfa,xtest,suppo,method)
  for (k in 1:tuldim){
     if (pc$frek[k]>0) apu[i,k]<-pc$resu[k]/pc$frek[k]
  }
  i<-i+1
}
# means std:s of apu:s rows
#tu<-meanstd(apu)
for (lo in 1:tuldim){
#   for (ek in 1:wv){
#      cv[lo]<-cv[lo]+apu[ek,lo]
#   }
#   cv[lo]<-cv[lo]/wv
#   for (ek in 1:wv){
#      cvstd[lo]<-cvstd[lo]+(apu[ek,lo]-cv[lo])^2
#   }
#   cvstd[lo]<-sqrt(cvstd[lo]/wv)
#   #
cv[lo]<-mean(apu[,lo],na.rm=T)              
cvstd[lo]<-sd(apu[,lo],na.rm=T)
for (vapp in 1:wv){
  if (!is.na(apu[vapp,lo])){
      numinmean[lo]<-numinmean[lo]+1
  }
}
}

return(list(alfa=pralfa,leafs=lnumber,cv=t(cv),cvstd=t(cvstd),numinmean=t(numinmean)))
#return(list(alfa=pralfa,leafs=lnumber,cv=t(tu$cv),cvstd=t(tu$cvstd)))
}






