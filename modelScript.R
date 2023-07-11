
# Multispecies occupancy model for two or more interacting species (Rota et al., 2016 Methods in Ecology and Evolution)
# Adapted from Kery and Royle, 2021 AHM 2
# Ugyen Penjor, October 2021 with generous support from Late Mike Meredith

# Load packages
library(AHMbook)
library(abind)
library(jagsUI)

# Load data
load("interactionModelData_Penjor2021_dcat.RData")
ls()

# Check data 
str(dataList1)

# Species detection history
str(ylist)

# Species detection as an array (site x rep x species)
str(y)

# Get sample sizes - constants
( nsites <- dim(y)[1] )			# sites
( nsurveys <- dim(y)[2] )		# surveys
( nspec <- dim(y)[3] )			# species
( ncat <- 2^nspec )				  # possible community states

# Prepare site covariates for occupancy (note that we are standardising using the scale() function)
settlement <- scale(dataList1$covars[, 'set'])
forest <- scale(dataList1$covars[, 'forest'])
sambar <- scale(dataList1$covars[, 'sam'])
gaur <- scale(dataList1$covars[, 'gau'])
serow <- scale(dataList1$covars[, 'ser'])

# Detection covariates
Disturb <- scale(dataList1$covars[, 'disturb'])
Trail <- dataList1$covars[, 'trail']

# Convert to matrix
Trail <- matrix(cbind(rep(Trail, 4)), ncol=4)
Disturb <- matrix(cbind(rep(Disturb, 4)), ncol=4)

# Condense multi-species detection array to be site x survey
ycat <- apply(y, c(1,2), paste, collapse="")
ycat[ycat == "000"] <- 1               # Unoccupied ('U')
ycat[ycat == "100"] <- 2               # Only tiger ('T') detected
ycat[ycat == "010"] <- 3               # Only leopard ('L') detected
ycat[ycat == "001"] <- 4               # Only dhole ('D') detected
ycat[ycat == "110"] <- 5               # Tiger and leopard ('TL')
ycat[ycat == "011"] <- 6               # Leopard and dhole ('LD')
ycat[ycat == "101"] <- 7               # Tiger and dhole ('TD')
ycat[ycat == "111"] <- 8               # All three species ('TLD')

# Convert each column to a numeric: this is our response variable
ycat <- apply(ycat, 2, as.numeric)     # check for warnings produced by NAs in observation

# Prepare covariate 
# Marginal occupancy covariates
psi_cov <- matrix(NA, ncol=6, nrow=nsites)
psi_cov[, 1] <- 1					# Intercept
psi_cov[, 2] <- settlement				
psi_cov[, 3] <- forest
psi_cov[, 4] <- sambar
psi_cov[, 5] <- gaur
psi_cov[, 6] <- serow

# Interaction covariates
psi_inxs_cov <- matrix(NA, ncol=2, nrow=nsites)
psi_inxs_cov[, 1] <- 1					# Intercept
psi_inxs_cov[, 2] <- settlement			

# Detection covariate
rho_cov <- array(NA, dim=c(nsites, nsurveys, 3))
rho_cov[,, 1] <- 1					# Intercept
rho_cov[,, 2] <- Trail	
rho_cov[,, 3] <- Disturb

# Detection - interaction
rho_inxs_cov <- rep(1, nsites)		# Intercept only

# Bundle and summarise data set
str( bdata <- list(y=ycat, psi_cov=psi_cov,
                   psi_inxs_cov=psi_inxs_cov, 
                   rho_cov=rho_cov, 
                   rho_inxs_cov=rho_inxs_cov, 
                   nsites=nsites, nsurveys=nsurveys, nfirst_order_psi=ncol(psi_cov), nsecond_order_psi=ncol(psi_inxs_cov), 
                   nfirst_order_rho=dim(rho_cov)[3], nsecond_order_rho=1, ncat=ncat) )


# Set initial values
# Get the maximum possible states across all potential surveys at a site
zinit <- apply(y, c(1,3), sum, na.rm=TRUE)
zinit[zinit > 1] <- 1          # convert to binary

# Convert to a category
zcat <- apply(zinit, 1, paste, collapse = "")
zcat[zcat == "000"] <- 1       # no species here
zcat[zcat == "100"] <- 2       # only tiger 
zcat[zcat == "010"] <- 3       # only leopard
zcat[zcat == "001"] <- 4       # only dhole
zcat[zcat == "110"] <- 5       # tiger and leopard
zcat[zcat == "011"] <- 6       # leopard and dhole
zcat[zcat == "101"] <- 7       # tiger and dhole
zcat[zcat == "111"] <- 8       # all three present

# Make numeric again
zcat <- as.numeric(zcat)

# Inits function
inits <- function() list(z=zcat)

# Parameters wanted
params <- c('betaT', 'betaL', 'betaD', 'betaTL', 'betaTD', 'betaLD', 'alphaT', 'alphaL', 'alphaD', 'mean.psiT', 'mean.psiL', 'mean.psiD', 'mean.pT', 'mean.pL', 'mean.pD', 'z')

# Call JAGS ~ 10 hrs for 1000 iter (!) - so you need a good machine to run longer iterations
outMod <- jags(bdata, inits, params, 'intModel.txt', n.chains=3, n.adapt=5000, n.burnin=1000, n.iter=50000, n.thin=20, parallel=TRUE)

print(outMod)

##############################################################################################################
