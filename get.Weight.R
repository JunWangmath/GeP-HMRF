###### the function get.Weight of each gene ######
###### JUN WANG, PKU ######

get.Weight=function(ppi.matrix.noself.6,sherlock.bf.5570.logLR,meth){
  if(meth=="eq1"){
    weight.gene.6=matrix(0.5,nrow=1,ncol=ncol(ppi.matrix.noself.6))
    print(meth)
  }else{
    if(meth=="sher"){
      print(meth)
      point1.x.weight=-2
      point1.y.weight=0.1
      point2.x.weight=1
      point2.y.weight=1
      print(paste0("uning sherlock log(BF) as weight mapping: ","(",point1.x.weight," ",point1.y.weight,"), (",point2.x.weight," ",point2.y.weight,")"))
      
      beta.in.weight=(point2.y.weight-point1.y.weight)/(point2.x.weight-point1.x.weight)
      alpha.in.weight=point1.y.weight - point1.x.weight*beta.in.weight
      
      weight.gene.6=sherlock.bf.5570.logLR*beta.in.weight+alpha.in.weight
      weight.gene.6[which(weight.gene.6<=point1.y.weight)]=point1.y.weight
      weight.gene.6[which(weight.gene.6>=point2.y.weight)]=point2.y.weight
      #write.csv(weight.gene.6,"weight.gene.6.csv",row.names=T)
    }else{
      print("weight by degree")
      weight.gene.6=apply(ppi.matrix.noself.6,2,sum)^0.5
      #write.csv(weight.gene.6,"weight.gene.6.csv",row.names=T)
    }
  }
  return(weight.gene.6)
}