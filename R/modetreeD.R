modetreeD<-function(d,treeseq,indseq,suppo,tot=1000){
#
alkm<-length(indseq)
xcoor<-matrix(0,tot,d)
ycoor<-matrix(0,tot,1)
zcoor<-matrix(0,tot,1)
mlabel<-matrix(0,tot,1)
#
laskuri<-1
for (i in alkm:1){
     inds<-indseq[i]
     leafnum<-treeseq$leafs[inds]
     alpha<-treeseq$alfa[inds]
     tree<-picktree(treeseq,leafnum)
     pv<-partition(tree,suppo)
     recs<-pv$recs
     values<-pv$values
     pg<-profgene(values,recs,frekv=F,cvol=T,ccen=T,cfre=F)
     mlkm<-moodilkm(pg$parent)
     for (j in 1:mlkm$lkm){  
        loca<-mlkm$modloc[j]
        for (dim in 1:d){ 
          xcoor[laskuri,dim]<-pg$center[dim,loca]
        }
        ycoor[laskuri]<-leafnum
        zcoor[laskuri]<-alpha
        mlabel[laskuri]<-j
        laskuri<-laskuri+1
     }
     mlkmpre<-mlkm
}
xcoor<-xcoor[1:(laskuri-1),]
ycoor<-ycoor[1:(laskuri-1)]
zcoor<-zcoor[1:(laskuri-1)]
mlabel<-mlabel[1:(laskuri-1)]
#plot(xcoor,ycoor,col=colo)
return(list(xcoor=xcoor,ycoor=t(ycoor),zcoor=t(zcoor),mlabel=t(mlabel)))
}                               








