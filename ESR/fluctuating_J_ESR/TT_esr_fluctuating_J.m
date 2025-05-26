%% Simulate quintet ESR spectrum using the fluctuating-J mechanism
clear
nrm = @(xs) xs/max(abs(xs));
load('ESR_1.5us.mat')% load data if you want

grd = sphgrid('Dinfh',61);% use Dinfh if coplanar, takes ages otherwise
weights = grd.weights;vecs = grd.vecs;
ngrd = length(weights);

% Set spin parameters. jhigh/jlow are the high and low exchange values
% k_flap is the rate of fluctuation b/w jlw/jhigh
% D,E are ZFS parameters
% bvec is B field vector; probably don't change this
jlow = 7.5e4;jhigh = 1.25e5;% 1e5 is 100GHz
k_flap = 3*(jhigh+jlow)/2;
D=1350;E=0;% Don't need to set Euler angles, that's later
bvec=[0,0,350];

% Set rates, same as above
k_tt = 0;
ks = [k_tt,k_flap];
ts = [1.5];% the timeslice we want, in microseconds

% need to define spin system to find the resfields
sys.S=[1 1];sys.D=[D,E;D,E];sys.J=(jlow+jhigh)/2;
eul1 = [0,0*pi/180,0];eul2 = [0,0*pi/180,0];
R1 = erot(eul1);R2 = erot(eul2);
sys.DFrame=[eul1;eul2];sys.lwpp=2;

Exp.mwFreq = 9.704081;
Exp.Range = [d(1,1),d(end,1)];
Exp.Harmonic = 0;
Exp.nPoints = length(d(:,1));

specs = zeros(1,Exp.nPoints);
Exp.MolFrame = [0,0,0];
% We loop over different molecule orientations, keeping B||z
f = waitbar(0,'Please wait...');
for k = 1:ngrd
    vec = vecs(:,k);
    [phi,theta] = vec2ang(vec);
    euls = [phi,theta,phi];
    eul1 = eulang(R1*erot(euls));
    eul2 = eulang(R2*erot(euls));
    sys.DFrame = [eul1;eul2];
    X = [D,E,D,E,jlow,jhigh,eul1,eul2,[0,0,350]];
    rho = fluctuating_J_populations(X,"pulse",ks,ts);%rho is TT density matrix
    pops = diag(rho);
    sys.initState = {pops,'eigen'};
    [B,spec] = pepper(sys,Exp);
    specs = specs + weights(k)*spec(1,:);
    waitbar(k/ngrd,f,"Loading EPR spectrum "+num2str(100*k/ngrd)+"%");
end
close(f)

% --- Begin optional bit to compare to just 1TT pop
clear Exp% get rid of molFrame
Exp.mwFreq = 9.7;
Exp.Range = [d(1,1),d(end,1)];
Exp.Harmonic = 0;
Exp.nPoints = length(d(:,1));

sys.initState = {[0,0,0,0,0,0,0,0.0,1],'coupled'};
[B,tt1spec] = pepper(sys,Exp);
% --- End optional
%% plot TAT ESR
figure('Position',[147.2222 131.2222 747.1111 450.2222])
plot(B,nrm(tt1spec),'LineWidth',2)
hold on
plot(d(:,1),nrm(d(:,2)),'k.','MarkerSize',12)
plot(B,nrm(specs(1,:))-1.2,'LineWidth',2)
plot(d(:,1),nrm(d(:,2))-1.2,'k.','MarkerSize',12)
axis tight
yticks([])
xlabel("Magnetic Field (mT)")
ylabel("Normalised ESR Spectrum")
leg1 = legend("Overlap model","","Fluctuating-J model","Data",'Location','west');
set(leg1,'Box','off')
fontsize(18,"points")
hA = gca;
hA.XRuler.MinorTickValues = [270:20:410];
hA.LineWidth = 2;hA.XMinorTick='on';
outerpos = hA.OuterPosition;
ti = hA.TightInset; 
left = outerpos(1) + ti(1);
bottom = outerpos(2) + ti(2);
ax_width = 0.99*(outerpos(3) - ti(1) - ti(3));
ax_height = outerpos(4) - ti(2) - ti(4);
hA.Position = [left bottom ax_width ax_height];



















