# Script created by Sam Hince
# 01/11/2021

##########################################################################

# nomenclature:
# V = velocity of aircraft
# phi = 
#
#
#
# B = number of blade
#
#
#

##########################################################################


## define constants 

V <- 100 # tf / s volocity

## calculate relevant starting conditions

## loop

# equation 1: maybe dont loop? 

phi <- atan((V*(1+a))/(Omega*R(1-aprime))) # should this be atan2? 

phi <- atan(V/(Omega*R)) # initial phi condition (step 1)

# quations 2:
	
alpha <- beta - phi

W <- V * (1 + a) * sin(phi)
RN <- rhp * W * (c / mu)

# equation 3: 

Cy <- (Cl * cos(phi)) - (Cd * sin(phi))

Cx <- (Cl * sin(phi)) + (Cd * cos(phi))

# equation 4:


