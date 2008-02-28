profavah<-function(dendat,
B,leaf,minlkm,seed,
sample="baggworpl",
prune="off",
splitscan=0,
seedf=1,
Q=NULL,
src="c",method="loglik",
frekv=NULL,cvol=TRUE,ccen=TRUE,cfre=FALSE)
{
#set.seed(seed=1)
#dendat<-matrix(rnorm(20),10)
#h<-1 
#N<-c(4,4)
#Q<-3

d<-dim(dendat)[2]

tr<-densfore(dendat,B,leaf,minlkm,seed,
sample,prune,splitscan,seedf,src=src,method=method)
left<-tr$left
right<-tr$right
val<-tr$val
vec<-tr$vec
nodenumOfTree<-length(left)

# make parent
parent<-matrix(0,length(left),1)
node<-1
while (node<=length(left)){
   if ((!is.na(left[node])) && (left[node]!=0)){
        parent[left[node]]<-node
   }
   if ((!is.na(right[node])) && (left[node]!=0)){
        parent[right[node]]<-node
   }
   node<-node+1
}

mi<-makeinfo(tr$left,tr$right,tr$mean,tr$low,tr$upp)
infopointer<-mi$infopointer
terminalnum<-mi$terminalnum
low<-mi$low
upp<-mi$upp
nodefinder<-mi$nodefinder
value<-mi$value

{
if (!is.null(Q)){
   maxval<-max(value)
   minval<-min(value)
   step<-(maxval-minval)/Q
   levseq<-seq(from=minval,to=maxval-step,by=step)
}
else{
   eppsi<-0.0000001
   levseq<-matrix(0,length(value),1)
   ordu<-order(value)
   ru<-1
   #car<-ordu[ru]
   #while ((ru <= length(value)) && (value[car]==0)){
   #     ru<-ru+1
   #     car<-ordu[ru]
   #}  # we have found first non zero
   laskuri<-1
   car<-ordu[ru]
   levseq[laskuri]<-value[car]-eppsi
   while (ru < length(value)){
       carnew<-ordu[ru+1]
       if (value[carnew]>value[car]){
          laskuri<-laskuri+1
          levseq[laskuri]<-value[carnew]-eppsi
      }
      ru<-ru+1
   }
   levseq<-levseq[1:laskuri]
   Q<-laskuri
}
}

levfrekv<-matrix(0,Q,1)
atomnum<-length(value)
for (i in 1:atomnum){
   for (j in 1:Q){
      if (value[i]>=levseq[j]){
         levfrekv[j]<-levfrekv[j]+1
      }
   }
}
numofall<-sum(levfrekv)

dentree<-decombag(numofall,levseq,
left,right,val,vec,infopointer,parent,
nodenumOfTree,terminalnum,
value,low,upp,nodefinder,
d)

invalue<-dentree$level
parent<-dentree$parent
component<-dentree$component
AtomlistAtom<-dentree$AtomlistAtom
AtomlistNext<-dentree$AtomlistNext

#if (cfre) nodefrek<-cfrekvdya(seplsets,binfrek) else nodefrek<-NULL 

if (ccen==TRUE) cvol<-TRUE
if (cvol){
  volume<-cvolumbag(component,AtomlistAtom,AtomlistNext,low,upp)
  kerroin<-cinte(invalue,volume,parent) 
  sepvalnor<-invalue/kerroin
} 
else{  
  volume<-NULL
  sepvalnor<-NULL
}

if (ccen && cvol){
  center<-ccentebag(component,AtomlistAtom,AtomlistNext,low,upp,volume)
  }
  else{
      center<-NULL
  }

return(list(parent=parent,level=sepvalnor,invalue=invalue,
volume=volume,center=center))    #nodefrek=nodefrek))
#values: normeeratut arvot
#invalues: alkuperaiset frekvenssit/arvot 
#nodefrek: kunkin solmun frekvenssi
}








