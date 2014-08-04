function output=main_file(input)

%main_file providing the calculations of the Fermi Problems assignment. Input required is: 
%input.n0 maximum density, input.R radius, input.RhoA (sets parameter A in the definition of rho(r)),
%input.RhoB (sets parameter B in the definition of rho(r)). Output provided is all cat. IV and II
%quantities: output.E0 minimum electric field max amplitude, output.I0 laser intensity, output.wl laser
%frequency, output.Vf focal volume, output.ncmax maximum packing factor, output.nfusjet number of fusion
%reactions in the cluster jet, output.ncluster number of clusters irradiated, output.vcluster volume of the
%cluster, output.rhoov contains parameters of the fitted rho profile, output.ClusterExpFus contains all
%parameters regarding cluster expansion and cluster fusion (see program GeneralCoulombModel_c)

%load cat. I quantities from the input structure:

R = input.R;
n0 = input.n0;
RhoA = input.RhoA;
RhoB = input.RhoB;

%initialisation cat.III quantities:

lambda=700e-9; %laser wavelength in metre cat
mi = 2*PhysConst.AMU; %see auxiliary file PhysConst.m
Z=1; %atom number 

rw = 1e-5; %laser waist in metre

%derive cat. IV quantities, which can be directly calculated

Vf = 2*((pi*(rw^2))^2)/(lambda); %focal volume
omegaL = (2*pi*PhysConst.c)/lambda; %laser frequency
Vcluster = (4/3) * pi * R^3; %volume cluster
Lray = pi*rw^2 / lambda ;

%initialse density profile
rvec=R/100:R/100:R; %create mesh
rhopiece = @(x,xdata) ((xdata<x(1)) + (xdata>=x(1)).*(xdata./x(1)).^(-x(2))).*(xdata<1);
rho = rhopiece([RhoA, RhoB],rvec./R);

%Number of particles per cluster
ncluster = 4 * pi * n0 * trapz(rvec, rvec.^2 .* rho);

%determine maximum packing factor, the function trapz(a,b) is the
%trapezoidal integral over \int_a(1)^a(end) b(a) da.
npack = ( omegaL^2 * PhysConst.epsilon0 * PhysConst.me * R^3) / (3*n0*PhysConst.e^2  * trapz(rvec, rvec.^2 .* rho));

%number of irradiated clusters
Nclusters = npack * (Vf / Vcluster);

%retrieve the size of the homogeneous part:
Rhom = RhoA.*R;

%Calculate the required electric field for electron ejection
omega0 = sqrt((PhysConst.e^2 * n0)/(3*PhysConst.epsilon0.*PhysConst.me)); %eigenfrequency oscillator
E0 = (2*Rhom*PhysConst.me)/PhysConst.e * abs(omegaL^2 - omega0^2); %mimimum required electric field
I0 = 0.5 * PhysConst.epsilon0 * PhysConst.c * E0^2 %electric field converted to intensity -> cat II unit, print it

%laser pulse energy, assuming a laser pulse duration of 20 fs - this
%parameter is required in order to check whether the result is physical or
%not
Ej=I0*pi*rw^2 * 20e-15 %cat II parameter -> print as output

%determine characteristic time of the expansion
t0 = sqrt(3.*mi.*PhysConst.epsilon0./((ncluster/Vcluster) .* Z * PhysConst.e^2));

%Calculate expansion and fusion reactions of the cluster:
expansion=ShockShellModelGeneralEuler(n0.*rho', R, Z, mi, 0:t0/20:10*t0, rvec', 1);

%Calculate total number of fusion reactions
Nfusion = expansion.FusY*Nclusters %cat II unit -> print

output.expansion=expansion;
output.Nfusion = Nfusion;
output.Ej = Ej;
output.I0 = I0;
output.E0 = E0;
output.Nclusters = Nclusters;
output.npack = npack;
output.ncluster = ncluster;
output.rho = rho;
output.rvec = rvec;
output.Vcluster = Vcluster;
output.Vf = Vf;
