%% This is the minimum working example for ESR of a triplet-pair, e.g in singlet fission and TTA
% written in MATLAB R2023b (should work in any version after R2021b) and Easyspin version 6.0.6
clear
% First we create the spin system called sys, and set values of various parameters
sys.S = [1,1];% two spin-1 particles, i.e a triplet-pair
% sys.D is a 2x2 matrix of ZFS paramters (MHz units)
% The top row is [D,E] for one triplet and bottom row is [D,E] for the other
% sys.D = [1414,14;1414,14] is roughly the TIPS-Tc parameters (see SI of 10.1103/PhysRevLett.125.097402) but don't sue me if it's wrong
sys.D = [1200,14;1200,14];
sys.J = 1e6;% (isotropic part of) inter-triplet exchange (MHz units)
sys.lwpp = 1.5;% Gaussian broadening parameter, Miles uses 1 to 2 (units of mT)
% sys.initState needs its own readme, sorry
sys.initState = {[0,0,0,0,1,0,1,0,1],'eigen'};

Exp.mwFreq = 9.95;
Exp.Range = [280,420];
Exp.Harmonic = 0;

Opt.separate = 'transitions';

[B,spec] = pepper(sys,Exp,Opt);
%%
figure('Position',[89 100 746 318])
plot(B,spec,LineWidth=2)
hold on
plot(B,sum(spec,1),'k',LineWidth=2)
xlabel('Field Strength (mT)')
ylabel('ESR Intensity (a.u.)')
