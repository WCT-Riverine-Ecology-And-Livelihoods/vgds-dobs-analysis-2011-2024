library(jagsUI)
library(rjags)
library(mcmcOutput)
library(dplyr)
library(ggplot2)

vgds_df <- read.csv("Data/vgdsDec2013_counts_rivcovs.csv")

vgds_df$straight_chan <- ifelse(vgds_df$channeltype == "S", 1, 0)
vgds_df$meander_chan <- ifelse(vgds_df$channeltype == "M", 1, 0)
vgds_df$island_chan <- ifelse(vgds_df$channeltype == "I", 1, 0)
vgds_df$log_rivwidth <- log10(vgds_df$rivwidth_mean)

for(i in 1:nrow(vgds_df)){
  vgds_df$total_y[i] <- sum(vgds_df$uniq_t1[i], vgds_df$uniq_t2[i], vgds_df$comm[i])
}

##Exploratory analysis
plot(vgds_df$total_y ~ vgds_df$W)
plot(vgds_df$total_y ~ vgds_df$G)
plot(vgds_df$total_y ~ vgds_df$F)
ggplot(data = vgds_df, aes(x = channeltype, y = total_y)) + geom_boxplot()

y <- as.matrix(vgds_df %>% select(uniq_t1, uniq_t2, comm))
colnames(y) <- NULL
data <- list(y = y, 
             M = nrow(vgds_df), #M - no. of sites
             W = vgds_df$W, 
             G = vgds_df$G,
             fog = vgds_df$F,
             rivwidth = vgds_df$log_rivwidth, 
             straight_chan = vgds_df$straight_chan,
             meander_chan = vgds_df$meander_chan,
             island_chan = vgds_df$island_chan)

cat("
model {
for(i in 1:M){ ## loop over M - no. of sites
# bionomial model for detection 
  logit(p[i,1]) <- alphaA0 + alpha2 * W[i] + alpha3 * G[i] +
                   alpha4 * rivwidth[i] + alpha5 * fog[i]
  logit(p[i,2]) <- alphaB0 + alpha2 * W[i] + alpha3 * G[i] +
                   alpha4 * rivwidth[i] + alpha5 * fog[i]

# poisson model for abundance 
  log(lambda[i]) <- beta0 + beta1 * straight_chan[i] + beta2 * meander_chan[i] +
                    beta3 * island_chan[i]

# probability of being detected by observer 1, observer 2, and both 
  pi[i,1] <- p[i,1] * (1-p[i,2]) 
  pi[i,2] <- p[i,2] * (1-p[i,1]) 
  pi[i,3] <- p[i,1] * p[i,2] 
  
# cumulative probability of being detected at all  
  pcap[i] <- pi[i,1] + pi[i,2] + pi[i,3]
  
# Poisson parameter:  multinomial cellprobs x expected abundance
  pi_lam[i,1] <- pi[i,1]*lambda[i]
  pi_lam[i,2] <- pi[i,2]*lambda[i]
  pi_lam[i,3] <- pi[i,3]*lambda[i]

# observed counts
  y[i,1] ~ dpois(pi_lam[i,1])
  y[i,2] ~ dpois(pi_lam[i,2])
  y[i,3] ~ dpois(pi_lam[i,3])
  
# missed animals  
  m[i] ~ dpois(lambda[i]*(1-pcap[i]))
  
# generate predictions of N_pred[i]  
  N_pred[i] ~ dpois(lambda[i])
  
# estimated realized abundance 
  N_realized[i] <- sum(y[i,1:3], m[i])
}
# Prior distributions (distribution - JAGS syntax)
# detection
pA0 ~ dunif(0,1) 
pB0 ~ dunif(0,1) 
alphaA0 <- logit(pA0) #intercept for obs1
alphaB0 <- logit(pB0) #intercept for obs2
# alpha1 ~ dnorm(0, 0.01) #coef for boatspeed
alpha2 ~ dnorm(0, 0.01) #coef for wind
alpha3 ~ dnorm(0, 0.01) #coef for glare
alpha4 ~ dnorm(0, 0.01) #coef for river width
alpha5 ~ dnorm(0, 0.01) #coef for fog

#abundance
beta0 ~ dnorm(0, 0.01)
beta1 ~ dnorm(0, 0.01) #for straight channel 
beta2 ~ dnorm(0, 0.01) #for meandering channel
beta3 ~ dnorm(0, 0.01) #for channel with islands

Npop <- sum(N_realized)
}
", fill = TRUE, file = "vgds_dobs_wo_boatspeed.txt")

# inits <- function(){
#   list (pA0 = runif(1), pB0 = runif(1),
#         alpha2 = runif(1), alpha3 = runif(1), alpha4 = runif(1),
#         beta0 = runif(1), beta1 = runif(1), beta2 = runif(1), beta3 = runif(1)
#         )
# }

# Define parameters to save and MCMC settings
params <- c("pA0", "pB0","alphaA0", "alphaB0", 
            "alpha2","alpha3", "alpha4", "alpha5",
            "beta0", "beta1", "beta2", "beta3", "Npop")
nc <- 3 ; ni <- 11000 ; nb <- 1000 ; nt <- 1

output <- jags(data, inits = NULL, params, 
               "vgds_dobs_wo_boatspeed.txt", 
               n.thin = nt, n.chains = nc, n.burnin = nb, n.iter = ni)

mco <- mcmcOutput(output)
diagPlot(mco)
plot(mco)
outputsum <- as.data.frame(summary(mco))

write.excel <- function(outputsum,row.names=TRUE,col.names=TRUE,...) {
  write.table(outputsum,"clipboard",sep="\t",row.names=row.names,col.names=col.names,...)
}

write.excel(outputsum)
