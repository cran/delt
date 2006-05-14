profaah<-function(dendat,B,leaf,minlkm=5,seed=1,
sample="bagg",prune="off",
splitscan=0,seedf=1,
scatter=0,
src="c",method="loglik")
{

hnum<-length(B)
hrun<-1
while (hrun<=hnum){
   hcur<-B[hrun]

   df<-
   densfore(dendat,hcur,leaf,minlkm=5,seed=1,
   sample="bagg",prune="off",
   splitscan=0,seedf=1,
   scatter=0,
   src="c",method="loglik")

   curtree<-proftree(df)

   if (hrun==1){
      if (hnum==1){
          treelist<-curtree
      }
      else{
          treelist=list(curtree)
      }
   }
   else{
      treelist=c(treelist,list(curtree))
   }
   hrun<-hrun+1
}
return(treelist)
}
