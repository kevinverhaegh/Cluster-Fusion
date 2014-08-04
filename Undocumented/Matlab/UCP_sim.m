clear all;
close all;

%%%%% filename parameters %%%%

paramdir = '100partLeffcs';
basedir = 'E:\GPT\NewSim\';  % all runs here

directosave = [basedir,'\',paramdir,'\'];
mkdir(directosave); 

e0=1.602e-19;
mu0=4*pi*10^-7;
auE=5.14220652e11;
aut=2.418884326505e-17;
auEn=4.35974417e-18;
auL=5.2917720859e-11;
auv=2.1876912633e6;

%Cluster parameters
n=1e27;					%Density (m^-3)
N=100;					%Amount of particles in a cluster
R0 = ((3*N)/(4*pi*n))^(1/3);	%Radius cluster
Ez= 15.5*PhysConst.e;			%Ionization energy (assume single ionized atoms)

%Laser parameters
tFWHM=50.0e-15;				%Full-Width Half Maximum laser

I0 = 5.0e21;				%Maximum laser-intensity
lambda = 800.0e-9;			%Wavelength laser
E0=sqrt((2*I0)/(PhysConst.c*PhysConst.epsilon0));	%Maximum electric field laser
omega=2*pi*PhysConst.c/lambda;		%Laser frequency
wd=omega;
tau = tFWHM/(2.0*sqrt(log(2.0)));	%Laser time scaling

t2=10*tau;
thalf = t2/2;

t1 = thalf;  % fine time step cutoff  (step 1)
step1 = 0.01*(1/omega);   % fine step size X divided by the plasma frequency
step2 = 0.1*(1/omega);   % coarse step size X divided by the plasma

tvec = [0:step1:(t1-step1) t1:step2:thalf];
env.A = E0;
env.omega=wd;
env.tau = tau;
env.type='gaussian';
Efield = wd*enveta(wd.*(tvec-thalf), env).*cos(wd.*(tvec-thalf));
for i=1:numel(tvec)
    if abs(Efield(i)/auE) > (((Ez/auEn)^2)/4)
        break
    end
end

tion=tvec(i);

IonModel = 0;                         %Ionization model, 0 - instantaneous, 1 - linear, 2 - Landau, 3 - ADK
%Assume center cluster z=0, ionization time, assume instantaneous
%ionization, no penetration depth
Field = 1;                          %0 - no field, 1 - plane wave gaussian, 2 - plane wave cos^2, 3-.., 5-Quesnel and Morra

IoniTim = 2*pi/wd;
if IonModel == 0
    IoniTim = 0;
end

if IonModel == 2
    LandauIonization(Field, E0, wd, tau, directosave, thalf,Ez)
    tion = 0;
    IoniTim = 0;
end

Eioncase = 0; %Energy distribution case for electrons directly after ionization 0 for no energy, 1 for ADK prediction distribution linear polarization,

whichmodel = 1;
SpaceCharge = 1;

%if this value is no integer, maxwell distribution is assumed with this
%value the temperature in eV

%%%%% parameters %%%%
mr.SpaceCharge = SpaceCharge; %0 for none, 1 for electron-electron, 2 for all
mr.whichmodel = whichmodel;  %0 for Spacecharge3DTree, 1 for Spacecharge3D, 2 for Spacecharge3DMesh model
mr.nps = N; % number of electrons/ions; total number of simulated particles equals 2*nps
mr.dens = n; % electron density in [m^-3]
mr.rseed = -1; % random seed; different rand. seeds cause different initial pos and velocities, rseed=-1 will cause the random seed to be the current time
mr.eps = 0.01; % Coulomb potential roud-off parameter
mr.wd = omega; % angular frequency of the external electric field
mr.E0 = E0; % maximum electric field
mr.tau = tau;  %scaling parameter time
mr.acc = 10; % accuracy parameter for the Runge-Kutta solver
mr.theta = 0.9; % accuracy parameter for the Barnes-Hut algorithm
mr.verbose = 1;% show debug output
mr.Eioncase = Eioncase; %Case ion distribution
mr.Field = Field; %Field type
mr.IonModel = IonModel; %Ionization model used
mr.tion = tion; %Ionization time
mr.IoniTim = IoniTim; %Time in which ionization occurs
mr.Ez = Ez; %Ionization energy
mr.thalf = t2/2; %thalf

snapshot = 0; %snapshot = 0 shows 'tout' output, snapshot = 1 shows snapshots

cut_off = 0; % remove particles which are too far from the center of the sphere
whichParam = 'E0'; % parameter which will be varied

varyArray = [1e12, 0.5e12, 0.1e12, 0.05e12, 5e12, 1e13, 5e13] ;
%varyArray = sort(varyArray);
if strcmp(whichParam,'E0')==1
    varyArray2 = varyArray.*0;
end
len = length(varyArray);
Nvaried = 0;  % initialize Nvaried to 0, but change if N is varied
for i=1:len
    % redefine the varied parameter, for the ith output
    if strcmp(whichParam,'dens')
       dens = varyArray(i);
    elseif strcmp(whichParam,'nps')
        nps = varyArray(i);
    elseif strcmp(whichParam,'rseed')
        rseed = varyArray(i);
    elseif strcmp(whichParam,'eps')
        eps = varyArray(i);
    elseif strcmp(whichParam,'wp')
        wp = varyArray(i);
    elseif strcmp(whichParam,'theta')
        theta = varyArray(i);
    elseif strcmp(whichParam,'tau')
        tau = varyArray(i);
    elseif strcmp(whichParam,'Eioncase')
        T = varyArray(i);
    elseif strcmp(whichParam,'acc')      
        acc = varyArray(i);      
    elseif strcmp(whichParam,'E0')      
        E0 = varyArray(i);
        tvec = [0:step1:(t1-step1) t1:step2:thalf];
        env.A = E0;
        env.omega=wd;
        env.tau = tau;
        env.type='gaussian';
        Efield = wd*enveta(wd.*(tvec-thalf), env).*cos(wd.*(tvec-thalf));
        for j=1:numel(tvec)
            if abs(Efield(j)/auE) > (((Ez/auEn)^2)/4)
            break
            end
        end
        varyArray2(i) = tvec(j);
    elseif strcmp(whichParam,'wd')      
        wd = varyArray(i);
    elseif strcmp(whichParam,'IoniTim')
        IoniTim = varyArray(i);   
    end
end



w=strcat('mr.',whichParam);

% % Select the parameter to be varied:
% 

for i = 1:length(varyArray)
   
    if strcmp(whichParam,'dens')
       mr.dens = varyArray(i);
    elseif strcmp(whichParam,'nps')
        mr.nps = varyArray(i);
    elseif strcmp(whichParam,'rseed')
        mr.rseed = varyArray(i);
    elseif strcmp(whichParam,'eps')
        mr.eps = varyArray(i);
    elseif strcmp(whichParam,'wp')
        mr.wp = varyArray(i);
    elseif strcmp(whichParam,'theta')
        mr.theta = varyArray(i);
    elseif strcmp(whichParam,'tau')
        mr.tau=varyArray(i);
    elseif strcmp(whichParam,'Eioncase')
        mr.T = varyArray(i);
    elseif strcmp(whichParam,'acc')      
        mr.acc = varyArray(i); 
    elseif strcmp(whichParam,'E0')      
        mr.E0 = varyArray(i);
        mr.tion = varyArray2(i);
    elseif strcmp(whichParam,'wd')      
        mr.wd = varyArray(i);
    elseif strcmp(whichParam,'IoniTim')
        mr.IoniTim = varyArray(i);   
    end
    
  
    writeUCPmrfile(directosave,mr,i);
    
end

writeUCPinfile(directosave,cut_off,snapshot,t1,t2,step1,step2,whichmodel, SpaceCharge, IonModel, Eioncase);
writeUCPbatfile(directosave,length(varyArray), IonModel);

save([directosave,'parameters.mat'],'mr','varyArray','t1','t2','step1','step2','whichParam');


