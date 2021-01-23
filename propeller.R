# Script created by Sam Hince
# 01/11/2021


library(varhandle)

##########################################################################

# nomenclature:
# V = velocity of aircraft
# phi = 
# Omega = 
#
#
# B = number of blade
# R = radius along blade
#
#

##########################################################################


##########################################################################

## read in data files
setwd("/home/sam/Documents/classGitRepos/MAE195")
geom <- read.csv(file = './propSpecs/CESSNA150PROP.csv', header = FALSE)
coef <- read.csv(file = './propSpecs/NACA4415_RN500K_NCRIT9.csv', header = TRUE)

goem_general <- t(geom[1:8,1:2])
colnames(goem_general) = goem_general[1,]
goem_general <- as.data.frame(goem_general[])
goem_general <- goem_general[-1,]
goem_general <- unfactor(goem_general)
	
goem_station <- as.data.frame(geom[10:nrow(geom),1:3])
goem_station <- unfactor(goem_station)
colnames(goem_station) = goem_station[1,]
goem_station <- goem_station[-1,]
row.names(goem_station) <- 1:nrow(goem_station)
goem_station$R <- as.numeric(goem_station$R)
goem_station$CHORD <- as.numeric(goem_station$CHORD)
goem_station$BETA <- as.numeric(goem_station$BETA)

##########################################################################

mu <- goem_general$KIN_VISC * goem_general$DENSITY # N s/m^2 dynamic viscosity of the fluid

rho <- goem_general$DENSITY # kg/m^3 density of air

B <- goem_general$BLADES # number of blades

R <- goem_general$DIAMETER / 2 # radius along the prop

V <- goem_general$VELOCITY # m / s volocity

Omega <- (2 * pi * goem_general$RPM) / 60 # rad / s angular velocity of the propeller


##########################################################################
## ass-umptions

beta <- 20 * (pi / 180) # blade angle is this given???

## calculate relevant starting conditions

##########################################################################
# equation 1:
eqn_phi <- function(){
	phi <- atan2((V*(1+a)), (Omega*r(1-aprime))) # should this be atan2? 
	return(phi)
}

eqn_phi_estimate <- function(V, Omega, r){
	phi <- atan2(V, (Omega*r)) # initial phi condition (step 1)
	return(phi)
}

# quations 2:
eqn_alpha <- function(beta, phi){
	alpha <- beta - phi
	return(alpha)
}

eqn_W <- function(V, a, phi){
	W <- (V * (1 + a)) / sin(phi)
	return(W)
}

eqn_Re <- function(rho, W, c, mu){
	Re <- rho * W * (c / mu)
	return(Re)
}

# equation 3: 
eqn_Cx <- function(){
	Cx <- (Cl * sin(phi)) + (Cd * cos(phi))
	return(Cx)
}

eqn_Cy <- function(){
	Cy <- (Cl * cos(phi)) - (Cd * sin(phi))
	return(Cy)
}

# equation 4:
eqn_a <- function(){
	a <- ((sigma/(4*varF))*(Cy/(sin(phi)^2))) / (1+(sigma/(4*varF)*(Cy/(sin(phi)^2))))
	return()
}

eqn_aprime <- function(){
	aprime <- ((sigma/(4*varF))*(Cx/(sin(phi)*cos(phi)))) / (1+(sigma/(4*varF)*(Cy/(sin(phi)*cos(phi)))))
	return()
}

##########################################################################

for (station in 1:nrow(goem_station)){
	#separate our values for this loop
	r <- goem_station$R[station]
	c <- goem_station$CHORD[station]
	beta <- goem_station$BETA[station] * (pi/180) # convery to rad
	
	sprintf("Blade station: %g", r)
	
	#step 1
	phi <- eqn_phi_estimate(V, Omega, r)
	
	#step 2
	alpha <- eqn_alpha(beta, phi)
	# calculate RN here... should use a, but we dont have it? 
	
	#step 3
	phi_degrees <- phi * (180/pi)
	Cl <- approx(x = coef$ALPHA, y = coef$CL, xout = phi_degrees, method="linear")$y
	Cd <- approx(x = coef$ALPHA, y = coef$CD, xout = phi_degrees, method="linear")$y
	
	#step 4
	
	
}

phi <- eqn_phi_estimate(V, Omega, r)


# run eqn 2

# get cl and cd -- digitize plots and figure that out 

# run eqn 3

# run eqn 4

# run the other part of eqn 1

# repeat with:

phi <- phiOld + (0.4 * (phiNew - phiOld)) # 0.4 is variable, recopmended value by liebeck 
# until it converges 

# 
phit <- atan2((r/R)*tan(phi)) # atan or atan2? 
f <- (B/2)*((1-(r/R))/((sin(phit))^2))
varF <- (2/pi)acos(exp(-f)) # pi could be something else


## loop

signma <- (B*c)/(2*pi*r) # local solidity # used in step 6 

