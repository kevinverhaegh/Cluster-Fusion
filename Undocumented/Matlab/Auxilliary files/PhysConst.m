% Defines some basic physical constants.
classdef PhysConst
    properties (Constant)
        k_b = 1.3806503e-23; % Boltzmann constant [J*K^-1]
        AMU = 1.66053873e-27; % Atomic mass unit [kg]
        e = 1.60217646e-19; % Elementary charge [coulombs]
        epsilon0 = 8.854187e-12; % Vacuum permittivity [A^2*s^4*kg^−1*m^−3]
        mu0 = 4*pi*10^(-7); % Vaccum permeability [N*A^-2]
        me = 9.10938189e-31; % Mass of electron [kg]
        h_planck = 6.62606877e-34; %  Planck constant [J*s]
        h_bar = 1.054571726e-34; %reduced Planck constant h/2Pi [J*s]
        eV = PhysConst.e; % Electron volt [J]
        c = 2.99792458e8; % speed of light in vacuum [m s^-1]
        auv = 2.187e6; %atomic unit to speed conversion
        auE=5.14220652e11; %atomic unit to electric field conversion
        aut=2.418884326505e-17; %atomic unit to time conversion
        auEn=4.35974417e-18; %atomic unit to energy conversion
        auL=5.2917720859e-11; %atomic unit to length conversion
    end;
end