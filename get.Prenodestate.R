###### the function get.Prenodestate ######
###### A function to initialized the genes' state (configuration of the network) ######
###### JUN WANG, PKU ######
get.Prenodestate=function(ppi.matrix.noself.6,gwas.sig.with.genename){
  gene.pre=unique(gwas.sig.with.genename)
  #length(gene.pre)
  gene.pre=as.vector(gene.pre)
  #str(gene.pre)
  gene.6=colnames(ppi.matrix.noself.6)
  #str(gene.6)
  gene.6=as.vector(gene.6)
  gene.pre.6=intersect(gene.pre,gene.6)
  #length(gene.pre.6)
  #write.table(gene.pre.6,"gene.pre.txt")
  
  node.state=matrix(0,nrow=ncol(ppi.matrix.noself.6),ncol=1)
  for (i in 1:length(gene.pre.6)){
  	tmp13=which(gene.6==gene.pre.6[i])
  	node.state[tmp13]=1
  }
  print(paste("The number of state 1 nodes is",sum(node.state==1)))
  print(paste("The number of state 0 nodes is",sum(node.state!=1)))
  return(node.state)
}