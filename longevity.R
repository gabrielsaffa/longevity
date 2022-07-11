###########################################################################################################################################################################################
### R code for: Longevity of the solitary silvery mole-rat Heliophobius argenteocinereus and its implications for the assumed link between sociality and longevity in African mole-rats ###


setwd("C:/R_folder/longevity")

library(rethinking)
library(phytools)
library(ape)


## read and re-code the data
long_data <- read.csv("Longevity_data-update-only-Bath_FINAL.csv", header=TRUE)
spp <- long_data$species # species index
SC <- ifelse(long_data$social_system=="Social",2,1) # index sociality
SC1 <- ifelse(SC==2,1,0) # create an indicator for sociality as well so we can model it
SS <- ifelse(long_data$sample_size=="Small",1,ifelse(long_data$sample_size=="Medium",2,3)) # index sample size categories
SS <- (SS - mean(SS)) / (sd(SS)*2) # standardize sample size
SS_ex <- log(long_data$sample_size_exact) # log-transform exact sample size
SS_ex <- (SS_ex - mean(SS_ex,na.rm=TRUE)) / (sd(SS_ex,na.rm=TRUE)*2) # standardize exact sample size
BM <- (long_data$body_mass - mean(long_data$body_mass)) / (sd(long_data$body_mass)*2) # standardize body mass
ML_rc <- (long_data$ML_rc - mean(long_data$ML_rc)) / (sd(long_data$ML_rc)*2) # standardize recorded longevity
ML_rs <- (long_data$ML_rs - mean(long_data$ML_rs)) / (sd(long_data$ML_rs)*2) # standardize longevity residual

# read the tree
long_tree <- read.nexus("long_tree.nex")
long_tree <- ls.consensus(long_tree) # create a consensus tree 
long_tree$tip.label[long_tree$tip.label=="Fukomys_mechowi"] <- "Fukomys_mechowii" # rename to F. mechowii
tree_trimmed <- keep.tip(long_tree, spp) # combine tips of the phylogeny with species in the data
R_OU <- corMartins(1,phy=tree_trimmed,form=~spp,fixed=FALSE) # adjust the tree according to the OU model but let the alpha be estimated rather than fixed
V <- vcv(R_OU) # convert to a variance-covariance matrix
V <- V[spp,spp] # order the species
R <- V/max(V) # convert to a correlation matrix
N <- length(spp) # number of species
species <- 1:N # species index


## create a data list
d <- list(
  ML_rc=ML_rc,
  ML_rs=ML_rs,
  SC=SC,
  SC1=SC1,
  SS=SS,
  BM=BM,
  R=R,
  N=N,
  species=species
)


## maximum recorded longevity (ML_rc)
## we'll model the data in the form of structural causal model (SCM)
## that is, each variable will be modeled as a function of those variables with a direct arrow pointing to them, as implied by our DAG (see Supplementaray Figure 2)
## this way, we'll explicitly model the entire causal system, as represented in our DAG
## with SCM, the data are modeled in a joint posterior distribution, from which we'll then simulate each intervention
## we'll first estimate the total effect of sociality
## that is, we need to control for body mass and phylogeny, according to our DAG, to "close" non-causal (i.e. back-door) paths
# write the SCM in a list
m <- alist(
  # BM
  BM ~ multi_normal("rep_vector(0,N)",S), # since body mass is standardized, we declare the mean as the vector of zeros
  matrix[N,N]:S <- sigma_BM*R, # and estimate the phylogenetic effect on the covariance among species
  sigma_BM ~ exponential(1),
  
  # ML
  ML_rc ~ multi_normal(mu,S),
  mu <- a_ML[SC] + bBM_ML[SC]*BM, # longevity is a function of body mass
  matrix[N,N]:S <- sigma_ML*R, # as well as phylogeny
  a_ML[SC] ~ normal(0,1),
  bBM_ML[SC] ~ normal(0,0.5),
  sigma_ML ~ exponential(1)
)

# fit the model with ulam
MLrc_te <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## estimate the direct effect of sociality
## to estimate the direct effect, we need to control for sample size as a mediating variable
## that is, sample size is on the causal path between sociality and longevity: SC -> SS -> ML
# write the SCM in a list
m <- alist(
  # BM
  BM ~ multi_normal("rep_vector(0,N)",S),
  matrix[N,N]:S <- sigma_BM*R,
  sigma_BM ~ exponential(1),
  
  ## SC1
  SC1 ~ bernoulli(p),
  logit(p) <- a_SC + bBM_SC*BM + s[species], # sociality is a function of body mass
  vector[species]:s ~ multi_normal("rep_vector(0,N)",S), # and phylogeny, which we model as deviations from average log-odds of sociality in the sample, i.e. intercept
  matrix[N,N]:S <- sigma_SC1*R,
  a_SC ~ normal(0,1),
  bBM_SC ~ normal(0,0.5),
  sigma_SC1 ~ exponential(1),
  
  # ML
  ML_rc ~ multi_normal(mu,S),
  mu <- a_ML[SC] + bBM_ML[SC]*BM + bSS_ML[SC]*SS, # we condition on sample size as well
  matrix[N,N]:S <- sigma_ML*R,
  a_ML[SC] ~ normal(0,1),
  bBM_ML[SC] ~ normal(0,0.5),
  bSS_ML[SC] ~ normal(0,0.5),
  sigma_ML ~ exponential(1),
  
  # SS
  SS ~ normal(mu_SS,sigma_SS),
  mu_SS <- a_SS + bSC_SS*SC1, # sample size is a function of sociality
  a_SS ~ normal(0,1),
  bSC_SS ~ normal(0,0.5),
  sigma_SS ~ exponential(1)
)

# fit the model with ulam
MLrc_de <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## maximum longevity residual (ML_rs)
## estimate the total effect of sociality
## because maximum longevity residual already accounts for body mass, we need to control only for phylogeny to get the total effect
# write the SCM in a list
m <- alist(
  # ML
  ML_rs ~ multi_normal(mu,S),
  mu <- a_ML[SC], 
  matrix[N,N]:S <- sigma_ML*R,
  a_ML[SC] ~ normal(0,1),
  sigma_ML ~ exponential(1)
)

# fit the model with ulam
MLrs_te <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## estimate the direct effect of sociality
## to estimate the direct effect, we need to control for sample size as a mediating variable as well
# write the model in a list
m <- alist(
  # BM
  BM ~ multi_normal("rep_vector(0,N)",S),
  matrix[N,N]:S <- sigma_BM*R,
  sigma_BM ~ exponential(1),
  
  ## SC1
  SC1 ~ bernoulli(p),
  logit(p) <- a_SC + bBM_SC*BM + s[species], 
  vector[species]:s ~ multi_normal("rep_vector(0,N)",S), 
  matrix[N,N]:S <- sigma_SC1*R,
  a_SC ~ normal(0,1),
  bBM_SC ~ normal(0,0.5),
  sigma_SC1 ~ exponential(1),
  
  # ML
  ML_rs ~ multi_normal(mu,S),
  mu <- a_ML[SC] + bSS_ML[SC]*SS, 
  matrix[N,N]:S <- sigma_ML*R,
  a_ML[SC] ~ normal(0,1),
  bSS_ML[SC] ~ normal(0,0.5),
  sigma_ML ~ exponential(1),
  
  # SS
  SS ~ normal(mu_SS,sigma_SS),
  mu_SS <- a_SS + bSC_SS*SC1, 
  a_SS ~ normal(0,1),
  bSC_SS ~ normal(0,0.5),
  sigma_SS ~ exponential(1)
)

# fit the model with ulam
MLrs_de <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## use exact sample sizes instead of sample size categories to estimate direct effect of sociality
## because some species lack information on exact sample size, we'll impute those missing values while fitting the model
d$SS_ex <- SS_ex # add exact sample size in the data list

## maximum recorded longevity (ML_rc)
# write the SCM in a list
m <- alist(
  # BM
  BM ~ multi_normal("rep_vector(0,N)",S),
  matrix[N,N]:S <- sigma_BM*R,
  sigma_BM ~ exponential(1),
  
  ## SC1
  SC1 ~ bernoulli(p),
  logit(p) <- a_SC + bBM_SC*BM + s[species], # sociality is a function of body mass
  vector[species]:s ~ multi_normal("rep_vector(0,N)",S), # and phylogeny, which we model as deviations from average log-odds of sociality in the sample, i.e. intercept
  matrix[N,N]:S <- sigma_SC1*R,
  a_SC ~ normal(0,1),
  bBM_SC ~ normal(0,0.5),
  sigma_SC1 ~ exponential(1),
  
  # ML
  ML_rc ~ multi_normal(mu,S),
  mu <- a_ML[SC] + bBM_ML[SC]*BM + bSS_ML[SC]*SS_ex, # we condition on sample size as well
  matrix[N,N]:S <- sigma_ML*R,
  a_ML[SC] ~ normal(0,1),
  bBM_ML[SC] ~ normal(0,0.5),
  bSS_ML[SC] ~ normal(0,0.5),
  sigma_ML ~ exponential(1),
  
  # SS_ex
  SS_ex ~ normal(mu_SS,sigma_SS),
  mu_SS <- a_SS + bSC_SS*SC1, # sample size is a function of sociality
  a_SS ~ normal(0,1),
  bSC_SS ~ normal(0,0.5),
  sigma_SS ~ exponential(1)
)

# fit the model with ulam
MLrc_de_ex <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## maximum longevity residual (ML_rs)
# write the model in a list
m <- alist(
  # BM
  BM ~ multi_normal("rep_vector(0,N)",S),
  matrix[N,N]:S <- sigma_BM*R,
  sigma_BM ~ exponential(1),
  
  ## SC1
  SC1 ~ bernoulli(p),
  logit(p) <- a_SC + bBM_SC*BM + s[species], 
  vector[species]:s ~ multi_normal("rep_vector(0,N)",S), 
  matrix[N,N]:S <- sigma_SC1*R,
  a_SC ~ normal(0,1),
  bBM_SC ~ normal(0,0.5),
  sigma_SC1 ~ exponential(1),
  
  # ML
  ML_rs ~ multi_normal(mu,S),
  mu <- a_ML[SC] + bSS_ML[SC]*SS_ex, 
  matrix[N,N]:S <- sigma_ML*R,
  a_ML[SC] ~ normal(0,1),
  bSS_ML[SC] ~ normal(0,0.5),
  sigma_ML ~ exponential(1),
  
  # SS_ex
  SS_ex ~ normal(mu_SS,sigma_SS),
  mu_SS <- a_SS + bSC_SS*SC1, 
  a_SS ~ normal(0,1),
  bSC_SS ~ normal(0,0.5),
  sigma_SS ~ exponential(1)
)

# fit the model with ulam
MLrs_de_ex <- ulam(m, cores=4, chains=4, iter=2000, data=d)


## extract posterior samples
MLrc_te_post <- extract.samples(MLrc_te)
MLrc_de_post <- extract.samples(MLrc_de)
MLrs_te_post <- extract.samples(MLrs_te)
MLrs_de_post <- extract.samples(MLrs_de)
MLrc_de_ex_post <- extract.samples(MLrc_de_ex)
MLrs_de_ex_post <- extract.samples(MLrs_de_ex)


## for each variable in the SCM, we'll simulate a matrix with columns for species and rows for simulated observations
## we'll then loop over the entire matrix of simulated observations when simulating subsequent variables
## i.e. those that are functions of other variables as implied by our DAG
## while also averaging over the posterior distribution of each parameter of a given model
## having done that, we'll then compute the causal effect of sociality as a contrast, i.e. difference, in longevity between solitary and social species
## conceptually, it is: "E(Yx - Yx*)" for total effect and "E(Yxz - Yx*z)" for direct effect, where Y is the outcome, x and x* are levels of an intervening variable, and z is a specified value of a mediator (see Pearl 2001)
## to do that, we'll simulate observations of longevity as if all species were solitary (x) and social (x*), respectively
n_obs <- 4000 # number of posterior samples as well as observations to be simulated
n_spp <- 9 # number of species in the sample / to be simulated

## maximum recorded longevity - total effect
## we first simulate body mass
BM <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  BM[i,] <- rmvnorm(1,rep(0,n_spp),diag(MLrc_te_post$sigma_ML[i],n_spp,n_spp)%*%R) # third argument just replaces 1s at the diagonal of the correlation matrix R with estimated variances
}

## simulate counterfactual distributions of longevity, averaging over the distribution of body mass for each species
MLrc_te_SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
MLrc_te_SC2 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    MLrc_te_SC1[i,] <- rmvnorm(1,rep(MLrc_te_post$a_ML[i,1] + MLrc_te_post$bBM_ML[i,1]*BM[i,j],n_spp),diag(MLrc_te_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
    MLrc_te_SC2[i,] <- rmvnorm(1,rep(MLrc_te_post$a_ML[i,2] + MLrc_te_post$bBM_ML[i,2]*BM[i,j],n_spp),diag(MLrc_te_post$sigma_ML[i],n_spp,n_spp)%*%R) # social
  }
}

## compute the contrast by subtracting solitary from social species
MLrc_te_cont <- MLrc_te_SC2 - MLrc_te_SC1
mean(MLrc_te_cont) # 0.49
PI(MLrc_te_cont) # -0.98, 1.97


## maximum recorded longevity - direct effect
## simulate body mass
BM <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  BM[i,] <- rmvnorm(1,rep(0,n_spp),diag(MLrc_de_post$sigma_ML[i],n_spp,n_spp)%*%R)
}

## simulate sociality
p <- matrix(NA,nrow=n_obs,ncol=n_spp)
SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    p[i,j] <- inv_logit(MLrc_de_post$a_SC[i] + MLrc_de_post$bBM_SC[i]*BM[i,j] + MLrc_de_post$s[i,j])
    SC1[i,j] <- rbern(1,p[i,j])
  }
}

## simulate sample size
SS <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    SS[i,j] <- rnorm(1,MLrc_de_post$a_SS[i] + MLrc_de_post$bSC_SS[i]*SC1[i,j],MLrc_de_post$sigma_SS[i])
  }
}

## simulate counterfactual distributions of longevity, averaging over the distributions of body mass and sample size for each species
MLrc_de_SC1 <- matrix(NA,ncol=n_spp,nrow=n_obs)
MLrc_de_SC2 <- matrix(NA,ncol=n_spp,nrow=n_obs)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    MLrc_de_SC1[i,] <- rmvnorm(1,rep(MLrc_de_post$a_ML[i,1] + MLrc_de_post$bBM_ML[i,1]*BM[i,j] + MLrc_de_post$bSS_ML[i,1]*SS[i,j],n_spp),diag(MLrc_de_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
    MLrc_de_SC2[i,] <- rmvnorm(1,rep(MLrc_de_post$a_ML[i,2] + MLrc_de_post$bBM_ML[i,1]*BM[i,j] + MLrc_de_post$bSS_ML[i,2]*SS[i,j],n_spp),diag(MLrc_de_post$sigma_ML[i],n_spp,n_spp)%*%R) # social 
  }
}

## compute the contrast by subtracting solitary from social species
MLrc_de_cont <- MLrc_de_SC2 - MLrc_de_SC1
mean(MLrc_de_cont) # 0.43
PI(MLrc_de_cont) # -0.95, 1.78


## maximum longevity residual - total effect
## simulate counterfactual distributions of longevity
MLrs_te_SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
MLrs_te_SC2 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  MLrs_te_SC1[i,] <- rmvnorm(1,rep(MLrs_te_post$a_ML[i,1],n_spp),diag(MLrs_te_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
  MLrs_te_SC2[i,] <- rmvnorm(1,rep(MLrs_te_post$a_ML[i,2],n_spp),diag(MLrs_te_post$sigma_ML[i],n_spp,n_spp)%*%R) # social 
}

## compute the contrast by subtracting solitary from social species
MLrs_te_cont <- MLrs_te_SC2 - MLrs_te_SC1
mean(MLrs_te_cont) # 0.48
PI(MLrs_te_cont) # -0.95, 1.88


## maximum longevity residual - direct effect
## simulate body mass
BM <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  BM[i,] <- rmvnorm(1,rep(0,n_spp),diag(MLrs_de_post$sigma_ML[i],n_spp,n_spp)%*%R)
}

## simulate sociality
p <- matrix(NA,nrow=n_obs,ncol=n_spp)
SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    p[i,j] <- inv_logit(MLrs_de_post$a_SC[i] + MLrs_de_post$bBM_SC[i]*BM[i,j] + MLrs_de_post$s[i,j])
    SC1[i,j] <- rbern(1,p[i,j])
  }
}

## simulate sample size
SS <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    SS[i,j] <- rnorm(1,MLrs_de_post$a_SS[i] + MLrs_de_post$bSC_SS[i]*SC1[i,j],MLrs_de_post$sigma_SS[i])
  }
}

## simulate counterfactual distributions of longevity
MLrs_de_SC1 <- matrix(NA,ncol=n_spp,nrow=n_obs)
MLrs_de_SC2 <- matrix(NA,ncol=n_spp,nrow=n_obs)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    MLrs_de_SC1[i,] <- rmvnorm(1,rep(MLrs_de_post$a_ML[i,1] + MLrs_de_post$bSS_ML[i,1]*SS[i,j],n_spp),diag(MLrs_de_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
    MLrs_de_SC2[i,] <- rmvnorm(1,rep(MLrs_de_post$a_ML[i,2] + MLrs_de_post$bSS_ML[i,2]*SS[i,j],n_spp),diag(MLrs_de_post$sigma_ML[i],n_spp,n_spp)%*%R) # social 
  }
}

## compute the contrast by subtracting solitary from social species
MLrs_de_cont <- MLrs_de_SC2 - MLrs_de_SC1
mean(MLrs_de_cont) # 0.43
PI(MLrs_de_cont) # -0.94, 1.80


## simulate counterfactuals for direct effects models with exact sample sizes
## maximum recorded longevity - direct effect
## simulate body mass
BM <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  BM[i,] <- rmvnorm(1,rep(0,n_spp),diag(MLrc_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R)
}

## simulate sociality
p <- matrix(NA,nrow=n_obs,ncol=n_spp)
SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    p[i,j] <- inv_logit(MLrc_de_ex_post$a_SC[i] + MLrc_de_ex_post$bBM_SC[i]*BM[i,j] + MLrc_de_ex_post$s[i,j])
    SC1[i,j] <- rbern(1,p[i,j])
  }
}

## simulate sample size
SS <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    SS[i,j] <- rnorm(1,MLrc_de_ex_post$a_SS[i] + MLrc_de_ex_post$bSC_SS[i]*SC1[i,j],MLrc_de_ex_post$sigma_SS[i])
  }
}

## simulate counterfactual distributions of longevity, averaging over the distributions of body mass and sample size for each species
MLrc_de_ex_SC1 <- matrix(NA,ncol=n_spp,nrow=n_obs)
MLrc_de_ex_SC2 <- matrix(NA,ncol=n_spp,nrow=n_obs)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    MLrc_de_ex_SC1[i,] <- rmvnorm(1,rep(MLrc_de_ex_post$a_ML[i,1] + MLrc_de_ex_post$bBM_ML[i,1]*BM[i,j] + MLrc_de_ex_post$bSS_ML[i,1]*SS[i,j],n_spp),diag(MLrc_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
    MLrc_de_ex_SC2[i,] <- rmvnorm(1,rep(MLrc_de_ex_post$a_ML[i,2] + MLrc_de_ex_post$bBM_ML[i,1]*BM[i,j] + MLrc_de_ex_post$bSS_ML[i,2]*SS[i,j],n_spp),diag(MLrc_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R) # social 
  }
}

## compute the contrast by subtracting solitary from social species
MLrc_de_ex_cont <- MLrc_de_ex_SC2 - MLrc_de_ex_SC1
mean(MLrc_de_ex_cont) # 0.38
PI(MLrc_de_ex_cont) # -1.02, 1.82


## maximum longevity residual - direct effect
## simulate body mass
BM <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(i in 1:n_obs) {
  BM[i,] <- rmvnorm(1,rep(0,n_spp),diag(MLrs_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R)
}

## simulate sociality
p <- matrix(NA,nrow=n_obs,ncol=n_spp)
SC1 <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    p[i,j] <- inv_logit(MLrs_de_ex_post$a_SC[i] + MLrs_de_ex_post$bBM_SC[i]*BM[i,j] + MLrs_de_ex_post$s[i,j])
    SC1[i,j] <- rbern(1,p[i,j])
  }
}

## simulate sample size
SS <- matrix(NA,nrow=n_obs,ncol=n_spp)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    SS[i,j] <- rnorm(1,MLrs_de_ex_post$a_SS[i] + MLrs_de_ex_post$bSC_SS[i]*SC1[i,j],MLrs_de_ex_post$sigma_SS[i])
  }
}

## simulate counterfactual distributions of longevity
MLrs_de_ex_SC1 <- matrix(NA,ncol=n_spp,nrow=n_obs)
MLrs_de_ex_SC2 <- matrix(NA,ncol=n_spp,nrow=n_obs)
for(j in 1:n_spp) {
  for(i in 1:n_obs) {
    MLrs_de_ex_SC1[i,] <- rmvnorm(1,rep(MLrs_de_ex_post$a_ML[i,1] + MLrs_de_ex_post$bSS_ML[i,1]*SS[i,j],n_spp),diag(MLrs_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R) # solitary
    MLrs_de_ex_SC2[i,] <- rmvnorm(1,rep(MLrs_de_ex_post$a_ML[i,2] + MLrs_de_ex_post$bSS_ML[i,2]*SS[i,j],n_spp),diag(MLrs_de_ex_post$sigma_ML[i],n_spp,n_spp)%*%R) # social 
  }
}

## compute the contrast by subtracting solitary from social species
MLrs_de_ex_cont <- MLrs_de_ex_SC2 - MLrs_de_ex_SC1
mean(MLrs_de_ex_cont) # 0.36
PI(MLrs_de_ex_cont) # -1.02, 1.79


## plot posterior distributions of each effect and their contrast
## compute densities
MLrc_te_SC1_d <- density(MLrc_te_SC1)
MLrc_te_SC2_d <- density(MLrc_te_SC2)
MLrc_te_cont_d <- density(MLrc_te_cont)
MLrc_de_SC1_d <- density(MLrc_de_SC1)
MLrc_de_SC2_d <- density(MLrc_de_SC2)
MLrc_de_cont_d <- density(MLrc_de_cont)
MLrs_te_SC1_d <- density(MLrs_te_SC1)
MLrs_te_SC2_d <- density(MLrs_te_SC2)
MLrs_te_cont_d <- density(MLrs_te_cont)
MLrs_de_SC1_d <- density(MLrs_de_SC1)
MLrs_de_SC2_d <- density(MLrs_de_SC2)
MLrs_de_cont_d <- density(MLrs_de_cont)
MLrc_de_ex_SC1_d <- density(MLrc_de_ex_SC1)
MLrc_de_ex_SC2_d <- density(MLrc_de_ex_SC2)
MLrc_de_ex_cont_d <- density(MLrc_de_ex_cont)
MLrs_de_ex_SC1_d <- density(MLrs_de_ex_SC1)
MLrs_de_ex_SC2_d <- density(MLrs_de_ex_SC2)
MLrs_de_ex_cont_d <- density(MLrs_de_ex_cont)

## plot computed distributions
{
  
  par(mfrow=c(3,2), mar=c(5.1,5.1,1.1,1.1))
  
  ## maximum recorded longevity - total effect
  plot(MLrc_te_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="total effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("A", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrc_te_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrc_te_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrc_te_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrc_te_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
  
  ## maximum longevity residual - total effect
  plot(MLrs_te_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="total effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("B", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrs_te_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrs_te_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrs_te_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrs_te_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
  
  ## maximum recorded longevity - direct effect
  plot(MLrc_de_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="direct effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("C", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrc_de_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrc_de_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrc_de_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrc_de_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
  
  ## maximum longevity residual - direct effect
  plot(MLrs_de_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="direct effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("D", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrs_de_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrs_de_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrs_de_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrs_de_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
  
  ## maximum recorded longevity - direct effect (exact sample size)
  plot(MLrc_de_ex_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="direct effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("E", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrc_de_ex_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrc_de_ex_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrc_de_ex_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrc_de_ex_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
  
  ## maximum longevity residual - direct effect (exact sample size)
  plot(MLrs_de_ex_SC1_d, xlim=c(-4.5,4.5), ylim=c(0,1), xlab="direct effect", ylab="density", main="", col="lightseagreen",lwd=3,cex.lab=1.75, yaxt="n", xaxt="n", axes=FALSE)
  axis(1, at=seq(-4,4, length.out=5), labels=c("-4","-2","0","2","4"), cex.axis=1.75)
  axis(2, at=seq(0,1, length.out=5), labels=c("0","0.25","0.5","0.75","1"), cex.axis=1.75)
  mtext("F", 3, line=-1, adj=0.05, cex=1.5)
  lines(MLrs_de_ex_SC2_d, col="deepskyblue4", lwd=3)
  polygon(MLrs_de_ex_SC2_d, col=col.alpha("deepskyblue4", alpha=0.55), border=NA)
  polygon(MLrs_de_ex_SC1_d, col=col.alpha("lightseagreen", alpha=0.55), border=NA)
  lines(MLrs_de_ex_cont_d, col="gray75", lwd=3, lty=2)
  
  # add legend
  legend <- c("solitary","social","contrast")
  cols <- c("lightseagreen","deepskyblue4","gray75")
  legend(0.5,1, col=cols, legend=legend, lwd=3, lty=c(1,1,2), cex=1.75, box.col=NA)
  
}

dev.off()


################################################################################
################################################################################