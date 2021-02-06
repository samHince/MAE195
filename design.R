# Script created by Sam Hince
# 02/04/2021

##########################################################################

#givens
D
B
V
RPM
Airfoil
Thrust or power

Omega <- RPM * (1/60) * (2*pi)
mu <- goem_general$KIN_VISC * goem_general$DENSITY # N s/m^2 dynamic viscosity of the fluid

##########################################################################

# monenclature 
# delta = non dimensional veloicty 
# Omega = angular velocity in radians per second 

##########################################################################


#step 1 # starting value for delta
delta <- 0  

#step 2 # determine F and phi at each station
for(stations){
	phit <- atan((V/(Omega*R)) * (1 + (delta/2))) # atan or atan2? 
	
	f <- (B/2)*((1-(Xi)) / ((sin(phit))^2))
	varF <- (2/pi) * atan((exp(2*f)-1)^(1/2)) # Ideal version: varF <- (2/pi) * acos(exp(-f)) 
	
	# useing r*tan(phi) = R*tan(phit) = const
	phi <- atan2(R*tan(phit),r) # phi <- atan2(tan(phit),Xi)
}

#step 3
Wc <- (4*pi*(V^2)*G*delta)/(B*Omega*Cl)
Rn <- (rho * Wc) / mu

#step 4 # fixxxxxxxx
alpha_degrees <- alpha * (180/pi)
Cl <- approx(x = coef$ALPHA, y = coef$CL, xout = alpha_degrees, method="linear")$y
Cd <- approx(x = coef$ALPHA, y = coef$CD, xout = alpha_degrees, method="linear")$y
epsilon

#step 5
X <- (Omega * r)/V

a <- (delta / 2) * (cos(phi)^2) * (1 - (epsilon*tan(phi)))
aprime <- (delta / (2*X)) * cos(phi) * sin(phi) * (1 + (epsilon / tan(phi))) 

#step 6
W <- (V * (1+a))/sin(phi)

#step 7
c <- # from Wc but how? 
beta <- alpha + phi

#step 8
I1
I2
J1
J2




