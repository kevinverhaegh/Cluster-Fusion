disp('Input')
rhon=1e30
rhoExp=1e25
fil = rhoExp/rhon  
En=100e3 *PhysConst.e
tauL = 20e-15
Rw = 1e-5
lambda = 700e-9
v=sqrt(En/(0.5*PhysConst.AMU*2))
sigmavv = 1e-15 * 1e-6 / v



disp('Results')
E = ((2*(rhon^(3/5))*PhysConst.e)/(3*PhysConst.epsilon0)) * (((15*PhysConst.epsilon0 *En)/(4*pi*PhysConst.e^2))^(1/5))
R = ((15*PhysConst.epsilon0*En)/(4*pi*PhysConst.e^2 * rhon^2))^(1/5)

EnL = 0.5 *PhysConst.c *PhysConst.epsilon0 * pi * tauL * Rw^2 * E^2 
EnL = 0.5 *PhysConst.c *PhysConst.epsilon0 * pi * tauL * Rw^2 * (((2*(rhon^(3/5))*PhysConst.e)/(3*PhysConst.epsilon0)) * (((15*PhysConst.epsilon0 *En)/(4*pi*PhysConst.e^2))^(1/5)))^2

L = pi*Rw^2 / lambda
tau = L/v
V = 2*pi^2 * Rw^4 / lambda

fil = lambda/(2*pi*L)

Nf = V *rhon^2 * fil^2 * sigmavv * L
Nf = (2*pi^3 * Rw^6  / lambda^2) * rhon^2 * fil^2 * (1e-15 * 1e-6 / v)

NfE = Nf * 14.1e6 * PhysConst.e

Q=NfE/EnL

Qtry = 4.8e6 * (sigmavv/(En^(2/5))) * 14.1e6 * PhysConst.e * (lambda^2 * rhon^(4/5))/tauL 

PredTajima = tau*3e9

Q2 = (2/3) *((lambda^2) * rhon^(4/5) / tauL)  * ((1e-15 * 1e-6 / v)/(En^(2/5)) ) * 14.1e6 * PhysConst.e / (2* pi^2 * 0.5 *PhysConst.c *PhysConst.epsilon0  * (((2*PhysConst.e)/(3*PhysConst.epsilon0))^2)  * (((15*PhysConst.epsilon0)/(4*pi*PhysConst.e^2))^(2/5)))