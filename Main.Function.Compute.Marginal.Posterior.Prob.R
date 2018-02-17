###### computing marginal posterior probability using MCMC (Gibbs Sampling) for each gene ######
###### JUN WANG, PEKING UNIVERSITY, 20180216 ######

source("get.Weight.R")
source("get.Prenodestate.R")
source("MCMC.gibbs.speedup.R")
source("get.posterior.prob.R")
require(graphics)

args <- commandArgs(TRUE) #users should give 6 parameters to run
alpha.args=as.numeric(args[1]) #alpha parameter in paper, alpha = -5 is suggested
#str(h.args)
tao1.args=as.numeric(args[2]) #beta1 - beta2 parameter in paper, beta1 - beta2 = 0.8 is suggested
#str(tao1.args)
tao_1.args=as.numeric(args[3]) #beta2 - beta3 parameter in paper, beta2 - beta3 = -0.001 is suggested
#str(tao_1.args)
step.mcmc=as.numeric(args[4]) #number of Gibbs Sampling steps, 4000 steps is suggested
#str(step.mcmc)
task.id=as.numeric(args[5]) # if task.id is set as 1, the program will not permutate the log-likelihood 
                            # score computed from Sherlock, this will output the original marginal posterior
                            # probability. If task.id != 1, program will do permutation, to get the randomized
                            # marginal posterior probability. We suggest to do 2000 permutation, i.e. set task.id from 1 to 2000 in each run.
disease.num=as.numeric(args[6]) #which disease you want to run, see "disease.all" below

disease.all=c("AMD_Fritsche_2015_GTExV6_LD85","Barrett_08_GTExV6_LD85","Cholestreol_HDL_ONE_Eur_GTExV6_LD85","Cholestreol_TC_ONE_Eur_GTExV6_LD85","Cronhs.EU_Liu_2015_GTExV6_LD85","Franke_Cronhs_meta_GTExV6_LD85","HDL_Willer_2013_GTExV6_LD85","LDL_Willer_2013_GTExV6_LD85","TotalCholesterol_Willer_2013_GTExV6_LD85")
disease=disease.all[disease.num]
print(disease)

weight.meth="eq1"

input.filedir="./source.file.prepared.for.marginal.4coloc/"
output.filedir1=paste0("./output/output_marginal_probability_disease_id_",disease.num,"/")

if(!file.exists("./output/")){
  dir.create("./output/")
}

if(!file.exists(output.filedir1)){
  dir.create(output.filedir1)
}

output.file1=paste0(output.filedir1,"marginal.posterior.probability.gibbs_",disease,"_a_",alpha.args,"_b1_",tao1.args,"_b2_",tao_1.args,"_step_",step.mcmc,"_taskid_",task.id,".txt")

ppi.matrix.noself=read.table(paste0(input.filedir,"ppi.for.",disease,".20160102.newforLD.txt"),header=T,row.names=1)
#str(ppi.matrix.noself)
gene=read.table(paste0(input.filedir,"gene.name.for.",disease,".20160102.newforLD.txt"))
#str(gene)
sherlock.bf.5570=read.table(paste0(input.filedir,"sherlock.logLR.for.",disease,".20160102.newforLD.txt"),header=T)
#str(sherlock.bf.5570)
sherlock.bf.5570.logLR=sherlock.bf.5570$logLR

sherlock.all=read.table(paste0(input.filedir,"logLR.gene.LD..eqtl.gwas.",disease,".20160101_2.txt"),header=T)
str(sherlock.all)
sherlock.all.logLR=sherlock.all$logLR

index.gene.in.sherlock5570=match(sherlock.bf.5570$gene.name,sherlock.all$gene.name)
sum(is.na(index.gene.in.sherlock5570))


if(task.id>1){
  print("permutate sherlock logBF")
  sherlock.all.logLR=sherlock.all.logLR[sample(length(sherlock.all.logLR),length(sherlock.all.logLR),replace = FALSE)]
  sherlock.bf.5570.logLR=sherlock.all.logLR[index.gene.in.sherlock5570]
}
sum(sherlock.bf.5570.logLR>3)
sherlock.sig.with.genename=sherlock.bf.5570$gene.name[which((length(sherlock.bf.5570$gene.name)-rank(sherlock.bf.5570.logLR))<50)]
length(sherlock.sig.with.genename)
gwas.sig.with.genename=sherlock.sig.with.genename

node.state=get.Prenodestate(ppi.matrix.noself,gwas.sig.with.genename)
meth=weight.meth
weight.gene=get.Weight(ppi.matrix.noself,sherlock.bf.5570.logLR,meth)

pram1=c(alpha.args,tao1.args,tao_1.args)
permu.or.not=1

ptm <- proc.time()
if(meth=="eq1"){
  res.gibbs=MCMC.gibbs.speedup.V1(pram1,step.mcmc,permu.or.not,node.state,sherlock.bf.5570.logLR)
}else{
  res.gibbs=MCMC.gibbs.speedup.V2(pram1,step.mcmc,permu.or.not,node.state,sherlock.bf.5570.logLR)
}

### add back the single node
prob.single.all=1-1/(1+exp(sherlock.all.logLR+pram1[1]))
res.gibbs.single=prob.single.all[-index.gene.in.sherlock5570]
res.gibbs.bind=c(res.gibbs,res.gibbs.single)
gene.in.gibbs.result=c(sherlock.all$gene.name[index.gene.in.sherlock5570],sherlock.all$gene.name[-index.gene.in.sherlock5570])
res.gibbs.bind.gene=cbind(gene.in.gibbs.result,res.gibbs.bind)
colnames(res.gibbs.bind.gene)=c("Gene","Margin.Posterior.Prob")

write.table(res.gibbs.bind.gene,output.file1,quote=F,col.names=T,row.names=F)
print("successfully done")
proc.time() - ptm
