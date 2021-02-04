# Script created by Sam Hince
# 01/11/2021


library(varhandle)

##########################################################################

# nomenclature:
# V = velocity of aircraft
# phi = 
# Omega = 
# B = number of blade
# r = radius along blade

##########################################################################

convergence_minimum <- 0.00001

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
eqn_phi <- function(V, a, Omega, r, aprime){
	phi <- atan2((V*(1+a)), (Omega*r*(1-aprime))) # should this be atan2? 
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
eqn_Cx <- function(Cl, Cd, phi){
	Cx <- (Cl * sin(phi)) + (Cd * cos(phi))
	return(Cx)
}

eqn_Cy <- function(Cl, Cd, phi){
	Cy <- (Cl * cos(phi)) - (Cd * sin(phi))
	return(Cy)
}

# equation 4:
eqn_a <- function(sigma, varF, Cy, phi){
	a <- ((sigma/(4*varF))*(Cy/(sin(phi)^2))) / (1+(sigma/(4*varF)*(Cy/(sin(phi)^2))))
	return(a)
}

eqn_aprime <- function(sigma, varF, Cx, phi){
	aprime <- ((sigma/(4*varF))*(Cx/(sin(phi)*cos(phi)))) / (1+(sigma/(4*varF)*(Cx/(sin(phi)*cos(phi)))))
	return(aprime)
}

##########################################################################

df <- data.frame()

for (station in 1:nrow(goem_station)){
	#separate our values for this loop
	r <- goem_station$R[station]
	Xi <- r/R
	c <- goem_station$CHORD[station]
	beta <- goem_station$BETA[station] * (pi/180) # convery to rad
	
	# correction for final blade station:
	if(Xi >= 1){
		Xi <- 1 - 1e-15
	}
	
	print(sprintf("Blade station: %g", r))
	
	#step 1
	# initial assumptions: 
	phi <- eqn_phi_estimate(V, Omega, r)
	a <- 0
	aprime <- 0
	
	convergent <- FALSE
	loops_to_converge <- 0
	while(convergent == FALSE){
		#step 2
		alpha <- eqn_alpha(beta, phi)
		W <- eqn_W(V, a, phi)
		Re <- eqn_Re(rho, W, c, mu)
		
		#step 3
		alpha_degrees <- alpha * (180/pi)
		Cl <- approx(x = coef$ALPHA, y = coef$CL, xout = alpha_degrees, method="linear")$y
		Cd <- approx(x = coef$ALPHA, y = coef$CD, xout = alpha_degrees, method="linear")$y
		
		#step 4
		Cx <- eqn_Cx(Cl, Cd, phi)
		Cy <- eqn_Cy(Cl, Cd, phi)
		
		#step 5
		phit <- atan((Xi)*tan(phi)) # atan or atan2? 
		f <- (B/2)*((1-(Xi)) / ((sin(phit))^2))
		varF <- (2/pi) * atan((exp(2*f)-1)^(1/2)) # Ideal version: varF <- (2/pi) * acos(exp(-f)) 
		
		#step 6
		sigma <- (B*c)/(2*pi*r) # local solidity
		a <- eqn_a(sigma, varF, Cy, phi)				
		aprime <- eqn_aprime(sigma, varF, Cx, phi)
		
		# corrections for finite tip chord:
		if((abs(a) > .7) || (abs(aprime) > .7)){
			aprime = 0.4
		}
		
		#step 7
		phi_new <- eqn_phi(V, a, Omega, r, aprime)
		
		#step 8
		percent_change <- (abs(phi - phi_new)/phi)
		
		if(percent_change < convergence_minimum){ 	# if the percent change has decreased enough 
			convergent <- TRUE  					# then assume we have reached convergence
			print(sprintf("Phi converged in %g loops", loops_to_converge))
		} else {									# otherwise increment phi
			phi <- phi + (0.4 * (phi_new - phi))    # 0.4 is variable, recopmended value by liebeck 
			loops_to_converge <- loops_to_converge + 1
		}
	}
	
	#record data
	df <- rbind(df, data.frame(station, r, c, beta * (180/pi), alpha_degrees, phi * (180/pi), Re / 1000000, a, aprime, f, varF, W, Cx, Cy, Cl/Cd))
	
}

#reformat data
colnames(df) <- c("station", "r", "c", "beta", "alpha", "phi", "Re", "a", "aprime", "f", "F", "W", "Cx", "Cy", "L/D")
print(df)

#step 9
dr <- goem_station$R[2]-goem_station$R[1]

n <- goem_general$RPM / 60 

#Ct: 
coef_term <- (1/(rho * (n^2) * ((2*R)^4)))

integral <- 0
for(station in 1:(nrow(goem_station)-1)){
	#separate our values for this loop
	c <- goem_station$CHORD[station]
	W <- df$W[station]
	Cy <- df$Cy[station]

	d <- ((1/2) * rho * (W^2) * B * c * Cy * dr)
	integral <- integral + d
}

thrust <- integral
Ct <-  integral * coef_term

print((sprintf("Ct: %g", Ct)))

#Cp: 
coef_term <- (1/(rho * (n^3) * ((2*R)^5)))

integral <- 0
for(station in 1:(nrow(goem_station)-1)){
	#separate our values for this loop
	r <- df$r[station]
	c <- goem_station$CHORD[station]
	W <- df$W[station]
	Cx <- df$Cx[station]
	
	d <- ((1/2) * rho * (W^2) * B * c * Cx * Omega * r * dr) # Omega missing in the notes
	integral <- integral + d
}

Cp <- integral * coef_term

print((sprintf("Cp: %g", Cp)))


