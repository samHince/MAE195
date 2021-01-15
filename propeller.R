# Script created by Sam Hince
# 01/11/2021

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


## define constants 

V <- 100 # m / s volocity
Omega <- (2 * pi * 2000) / 60 # rad / s angular velocity of the propeller
R <- 1 # radius along the prop
B <- 2 # number of blades

# https://www.engineeringtoolbox.com/air-absolute-kinematic-viscosity-d_601.html
mu <- 18 * 10^(-6) # N s/m^2 dynamic viscosity of the fluid

rho <- 1.225 # kg/m^3 density of air

signma <- (B*c)/(2*pi*r)# local solidity # used in step 6 

## ass-umptions

beta <- 20 * (pi / 180) # blade angle is this given???

## calculate relevant starting conditions

# run eqn 1

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

# equation 1: maybe dont loop? 

phi <- atan2((V*(1+a))/(Omega*r(1-aprime))) # should this be atan2? 

phi <- atan2(V/(Omega*r)) # initial phi condition (step 1)

# quations 2:
	
alpha <- beta - phi

W <- V * (1 + a) * sin(phi)
Re <- rho * W * (c / mu)

# equation 3: 

Cy <- (Cl * cos(phi)) - (Cd * sin(phi))
Cx <- (Cl * sin(phi)) + (Cd * cos(phi))

# equation 4:

a <- ((sigma/(4*varF))*(Cy/(sin(phi)^2))) / (1+(sigma/(4*varF)*(Cy/(sin(phi)^2))))
aprime <- ((sigma/(4*varF))*(Cx/(sin(phi)*cos(phi)))) / (1+(sigma/(4*varF)*(Cy/(sin(phi)*cos(phi)))))
