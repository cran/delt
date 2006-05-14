dentseq<-function(dendat,minlkm){
#Makes a sequence of density trees
#
#dendat is n*d-matrix
#minlkm >=1 integer, we split so long as there is at least minlkm
#  observations in rectangle
#
#Returns
#
#dendat<-matrix(rnorm(20),10) on 10*2 matriisi
#
kanta<-supp(dendat,epsi=0) #support for density estimate
#
# Muodostetaan puitten jono
#
bt<-densplit(dendat,minlkm,kanta,blokki=50)     #bigtree
ini<-initial(bt$ssr,bt$left,bt$right)
bt$left<-ini$left
bt$right<-ini$right
treeseq<-pruseq(bt)
return(treeseq)
}



