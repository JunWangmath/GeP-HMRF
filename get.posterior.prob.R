##### when we get the configuration of a network, we need to compute the posterior ######
##### probability of the configuration, and then I will choose the largest one ######

##### 2014/05/27 ######
##### We think the posterior probability can be approximated by the log pseudo-likelyhood ######

get.posterior.prob=function(mcmc.result){
  pseudo.like=0
  node.state.tmp=mcmc.result[,1]
  log.odds=mcmc.result[,2]
  for(node.index in 1:length(node.state.tmp)){
  	if(node.state.tmp[node.index]!=1){
  	  pseudo.like=sum(pseudo.like,log((1+exp(log.odds[node.index]))))
  	}else{
  	  pseudo.like=sum(pseudo.like,log((1+exp(-log.odds[node.index]))))
  	}
  }
  pseudo.like=-pseudo.like
  return(pseudo.like)
}
