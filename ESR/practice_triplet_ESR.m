%% 
clear
sys.S = 1;
sys.D = [1414,14];
sys.lwpp = 1;
sys.initState = {[0,1,0],'eigen'};

Exp.mwFreq = 9.8;
Exp.Range = [260,420];
Exp.Harmonic = 0;

Opt.separate = 'transitions';

[B,spec] = pepper(sys,Exp,Opt);
%
%Plot the ESR spectrum
%
figure('Position',[89 500 746 318])
plot(B,spec,LineWidth=2)
xlabel('Field Strength (mT)')
ylabel('ESR Intensity (a.u.)')
