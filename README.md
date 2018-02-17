# GeP-HMRF V1.0

# Inferring Gene-Disease Association by an Integrative Analysis of eQTL GWAS and Protein-Protein Interaction data
#Jun WANG, Peking Univrsity, junwangmath@gmail.com

All the GWAS, eQTL and Protein-Protein Interaction dataset are stored at the "source.file.prepared.for.marginal.4coloc" folder. Unzip it and put it in the same folder of other R scripts. Unzip it and put it in the same folder of other R scripts.

The main function is "Main.Function.Compute.Marginal.Posterior.Prob.R", which has 6 parameters.

To run the main function:
   Rscript Main.Function.Compute.Marginal.Posterior.Prob.R -5 0.8 -0.001 4000 1 1
   
Six parameters which can be tune: ############################################################################################################################### alpha.args=as.numeric(args[1]) #alpha parameter in paper, alpha = -5 is suggested #str(h.args) tao1.args=as.numeric(args[2]) #beta1 - beta2 parameter in paper, beta1 - beta2 = 0.8 is suggested #str(tao1.args) tao_1.args=as.numeric(args[3]) #beta2 - beta3 parameter in paper, beta2 - beta3 = -0.001 is suggested #str(tao_1.args) step.mcmc=as.numeric(args[4]) #number of Gibbs Sampling steps, 4000 steps is suggested #str(step.mcmc) task.id=as.numeric(args[5]) # if task.id is set as 1, the program will not permutate the log-likelihood # score computed from Sherlock, this will output the original marginal posterior # probability. If task.id != 1, program will do permutation, to get the randomized # marginal posterior probability. We suggest to do 2000 permutation, i.e. set task.id from 1 to 2000 in each run. disease.num=as.numeric(args[6]) #which disease you want to run, see "disease.all" below

disease.all=c("AMD_Fritsche_2015_GTExV6_LD85","Barrett_08_GTExV6_LD85","Cholestreol_HDL_ONE_Eur_GTExV6_LD85","Cholestreol_TC_ONE_Eur_GTExV6_LD85","Cronhs.EU_Liu_2015_GTExV6_LD85","Franke_Cronhs_meta_GTExV6_LD85","HDL_Willer_2013_GTExV6_LD85","LDL_Willer_2013_GTExV6_LD85","TotalCholesterol_Willer_2013_GTExV6_LD85") #################################################################################################################################
