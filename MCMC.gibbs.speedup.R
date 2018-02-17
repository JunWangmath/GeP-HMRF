###### We want to use a MCMC version instead of the ICM version, the MCMC version will get the posterior distribution of gene's state Z by sampling ######
###### JUN WANG, PEIKING UNIVERSITY, 20151022 ##################################################################################################

###### The Vertion 2: using weight (non-equal weight),20151202

MCMC.gibbs.speedup.V2=function(pram,step.mcmc,permutate.order,node.state,sherlock.logLR){
  #node.state.track=matrix(100,nrow=length(node.state),ncol=step.mcmc)
  node.state.tmp=node.state
  row.sum=matrix(0,nrow=length(node.state.tmp),ncol=1)
  #node.neighbor.info=matrix(100000,nrow=length(node.state),ncol=2)
  pram.tmp=pram
  
  h.tmp=pram.tmp[1]
  tao.one.tmp=pram.tmp[2]
  tao.minus.one.tmp=pram.tmp[3]
  
  iteration.times=0
  
  if(permutate.order==0){
    while(iteration.times<step.mcmc){
      #print(iteration.times)
      #node.state.previous=node.state.tmp
      #log.odds.tmp=matrix(0,nrow=length(node.state.tmp),ncol=1)
      for (node.index in 1:length(node.state.tmp)){
        index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        betaij.one=tao.one.tmp*sum(weight.gene[index.state.one]+weight.gene[node.index])
        betaij.minus.one=tao.minus.one.tmp*sum(weight.gene[index.state.minus.one]+weight.gene[node.index])
        #betaij.one=tao.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #betaij.minus.one=tao.minus.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))        
        
        #Log.Odds.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        log.odds.node.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        #log.odds.node.i=Log.Odds(node.index,node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp,sherlock.logLR)
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #node.neighbor.info[node.index,1]=length(index.state.one)
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #node.neighbor.info[node.index,2]=length(index.state.minus.one)
        #log.odds.tmp[node.index,1]=log.odds.node.i
        #print(paste(paste("The log.odds of node",node.index),log.odds.node.i))
        #if(log.odds.node.i>0){
        #  node.state.tmp[node.index]=1
        #}else{
        #  node.state.tmp[node.index]=-1
        #}
        prob.1=1-1/(1+exp(log.odds.node.i))
        node.state.tmp[node.index]=rbinom(1,1,prob.1)
      }
      #p.change=Log.Odds.only.prior.speedup(node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp)
      #str(p.change)
      #rand.num=runif(length(node.state.tmp),0,1)
      #index.1=which(p.change>rand.num)
      #node.state.tmp[index.1]=1
      #node.state.tmp[-index.1]=-1
      #tmp.error=sum(abs(node.state.previous-node.state.tmp))/length(node.state.tmp)/2
      iteration.times=iteration.times+1
      row.sum=row.sum+node.state.tmp
      #print(paste("The iteration times is ",iteration.times))
      #print(paste("The temporary error now is ", tmp.error))
      #print(paste("The number of state 1 gene is ",sum(node.state.tmp==1)))
      #node.state.track[,iteration.times]=node.state.tmp
      #node.state.track[,4*iteration.times-2]=log.odds.tmp[,1]
      #node.state.track[,c(4*iteration.times-1,4*iteration.times)]=node.neighbor.info
      #node.info=cbind(node.state.tmp,log.odds.tmp)
      #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
    }
  }
  if(permutate.order==1){
    while(iteration.times<step.mcmc){
      permutation.order=sample(length(node.state),length(node.state))
      #node.state.tmp.previous=node.state.tmp
      #log.odds.tmp=matrix(0,nrow=length(node.state.tmp),ncol=1)
      for (index.pmt in 1:length(node.state.tmp)){
        node.index=permutation.order[index.pmt]
        index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        betaij.one=tao.one.tmp*sum(weight.gene[index.state.one]+weight.gene[node.index])
        betaij.minus.one=tao.minus.one.tmp*sum(weight.gene[index.state.minus.one]+weight.gene[node.index])
        #betaij.one=tao.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #betaij.minus.one=tao.minus.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        
        log.odds.node.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        #log.odds.node.i=Log.Odds(node.index,node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp,sherlock.logLR)
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #node.neighbor.info[node.index,1]=length(index.state.one)
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #node.neighbor.info[node.index,2]=length(index.state.minus.one)
        #log.odds.tmp[node.index,1]=log.odds.node.i
        #if(log.odds.node.i>0){
        #  node.state.tmp[node.index]=1
        #}else{
        #  node.state.tmp[node.index]=-1
        #}
        prob.1=1-1/(1+exp(log.odds.node.i))
        node.state.tmp[node.index]=rbinom(1,1,prob.1)
      }
      #p.change=Log.Odds.only.prior.speedup(node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp)
      #str(p.change)
      #rand.num=runif(length(node.state.tmp),0,1)
      #index.1=which(p.change>rand.num)
      #node.state.tmp[index.1]=1
      #node.state.tmp[-index.1]=-1
      #tmp.error=sum(abs(node.state.tmp.previous-node.state.tmp))/length(node.state.tmp)/2
      iteration.times=iteration.times+1
      row.sum=row.sum+node.state.tmp
      #print(paste("The iteration times is ",iteration.times))
      #print(paste("The error now is ", tmp.error))
      #print(paste("The number of state 1 gene is ",sum(node.state.tmp==1)))
      #node.state.track[,iteration.times]=node.state.tmp
      #node.state.track[,4*iteration.times-2]=log.odds.tmp[,1]
      #node.state.track[,c(4*iteration.times-1,4*iteration.times)]=node.neighbor.info
      #node.info=cbind(node.state.tmp,log.odds.tmp)
      #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
    }
  }
  #node.info=cbind(node.state.tmp,log.odds.tmp)
  #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
  #return(cbind(node.info,node.neighbor.info))
  #return(node.info)
  #row.ave=row.sum/step.mcmc
  return(row.sum/step.mcmc)
}

##########################################################################################################################################

###### The Vertion 1: using equal weight, 20151022
MCMC.gibbs.speedup.V1=function(pram,step.mcmc,permutate.order,node.state,sherlock.logLR){
  #node.state.track=matrix(100,nrow=length(node.state),ncol=step.mcmc)
  node.state.tmp=node.state
  row.sum=matrix(0,nrow=length(node.state.tmp),ncol=1)
  #node.neighbor.info=matrix(100000,nrow=length(node.state),ncol=2)
  pram.tmp=pram
  
  h.tmp=pram.tmp[1]
  tao.one.tmp=pram.tmp[2]
  tao.minus.one.tmp=pram.tmp[3]
  
  iteration.times=0
  
  if(permutate.order==0){
    while(iteration.times<step.mcmc){
      #print(iteration.times)
      #node.state.previous=node.state.tmp
      #log.odds.tmp=matrix(0,nrow=length(node.state.tmp),ncol=1)
      for (node.index in 1:length(node.state.tmp)){
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #betaij.one=tao.one.tmp*sum(weight.gene[index.state.one]+weight.gene[node.index])/2
        #betaij.minus.one=tao.minus.one.tmp*sum(weight.gene[index.state.minus.one]+weight.gene[node.index])/2
        betaij.one=tao.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        betaij.minus.one=tao.minus.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))        
        
        #Log.Odds.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        log.odds.node.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        #log.odds.node.i=Log.Odds(node.index,node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp,sherlock.logLR)
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #node.neighbor.info[node.index,1]=length(index.state.one)
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #node.neighbor.info[node.index,2]=length(index.state.minus.one)
        #log.odds.tmp[node.index,1]=log.odds.node.i
        #print(paste(paste("The log.odds of node",node.index),log.odds.node.i))
        #if(log.odds.node.i>0){
        #  node.state.tmp[node.index]=1
        #}else{
        #  node.state.tmp[node.index]=-1
        #}
        prob.1=1-1/(1+exp(log.odds.node.i))
        node.state.tmp[node.index]=rbinom(1,1,prob.1)
      }
      #p.change=Log.Odds.only.prior.speedup(node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp)
      #str(p.change)
      #rand.num=runif(length(node.state.tmp),0,1)
      #index.1=which(p.change>rand.num)
      #node.state.tmp[index.1]=1
      #node.state.tmp[-index.1]=-1
      #tmp.error=sum(abs(node.state.previous-node.state.tmp))/length(node.state.tmp)/2
      iteration.times=iteration.times+1
      row.sum=row.sum+node.state.tmp
      #print(paste("The iteration times is ",iteration.times))
      #print(paste("The temporary error now is ", tmp.error))
      #print(paste("The number of state 1 gene is ",sum(node.state.tmp==1)))
      #node.state.track[,iteration.times]=node.state.tmp
      #node.state.track[,4*iteration.times-2]=log.odds.tmp[,1]
      #node.state.track[,c(4*iteration.times-1,4*iteration.times)]=node.neighbor.info
      #node.info=cbind(node.state.tmp,log.odds.tmp)
      #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
    }
  }
  if(permutate.order==1){
    while(iteration.times<step.mcmc){
      permutation.order=sample(length(node.state),length(node.state))
      #node.state.tmp.previous=node.state.tmp
      #log.odds.tmp=matrix(0,nrow=length(node.state.tmp),ncol=1)
      for (index.pmt in 1:length(node.state.tmp)){
        node.index=permutation.order[index.pmt]
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #betaij.one=tao.one.tmp*sum(weight.gene[index.state.one]+weight.gene[node.index])/2
        #betaij.minus.one=tao.minus.one.tmp*sum(weight.gene[index.state.minus.one]+weight.gene[node.index])/2
        betaij.one=tao.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        betaij.minus.one=tao.minus.one.tmp*sum((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        
        log.odds.node.i=h.tmp+betaij.one+betaij.minus.one+sherlock.logLR[node.index]
        #log.odds.node.i=Log.Odds(node.index,node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp,sherlock.logLR)
        #index.state.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp==1))
        #node.neighbor.info[node.index,1]=length(index.state.one)
        #index.state.minus.one=which((ppi.matrix.noself[,node.index]!=0)&(node.state.tmp!=1))
        #node.neighbor.info[node.index,2]=length(index.state.minus.one)
        #log.odds.tmp[node.index,1]=log.odds.node.i
        #if(log.odds.node.i>0){
        #  node.state.tmp[node.index]=1
        #}else{
        #  node.state.tmp[node.index]=-1
        #}
        prob.1=1-1/(1+exp(log.odds.node.i))
        node.state.tmp[node.index]=rbinom(1,1,prob.1)
      }
      #p.change=Log.Odds.only.prior.speedup(node.state.tmp,h.tmp,tao.one.tmp,tao.minus.one.tmp)
      #str(p.change)
      #rand.num=runif(length(node.state.tmp),0,1)
      #index.1=which(p.change>rand.num)
      #node.state.tmp[index.1]=1
      #node.state.tmp[-index.1]=-1
      #tmp.error=sum(abs(node.state.tmp.previous-node.state.tmp))/length(node.state.tmp)/2
      iteration.times=iteration.times+1
      row.sum=row.sum+node.state.tmp
      #print(paste("The iteration times is ",iteration.times))
      #print(paste("The error now is ", tmp.error))
      #print(paste("The number of state 1 gene is ",sum(node.state.tmp==1)))
      #node.state.track[,iteration.times]=node.state.tmp
      #node.state.track[,4*iteration.times-2]=log.odds.tmp[,1]
      #node.state.track[,c(4*iteration.times-1,4*iteration.times)]=node.neighbor.info
      #node.info=cbind(node.state.tmp,log.odds.tmp)
      #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
    }
  }
  #node.info=cbind(node.state.tmp,log.odds.tmp)
  #print(paste("the pseodo.likelyhood is",get.posterior.prob(node.info)))
  #return(cbind(node.info,node.neighbor.info))
  #return(node.info)
  #row.ave=row.sum/step.mcmc
  return(row.sum/step.mcmc)
}
