clear
% initialise the triplet-pair however you please
sys.S = [1 1];
sys.D = [1000,100;1200,200];
eul1 = rand(1,3); eul2 = rand(1,3);
sys.DFrame = [eul1;eul2];
sys.J = 0;

% choose polar angles of field vector
phi = 2*pi*rand;theta = pi*rand;

figure
% first plot: use phi and theta to define field direction
subplot(2,1,1)
levelsplot(sys,[phi,theta],[-10,10])
hold on
title("These should be the same")
% second plot: keep B||z, and rotate molecule instead
R = erot([phi,theta,phi]);
sys.DFrame = [eulang(erot(eul1)*R');eulang(erot(eul2)*R')];
subplot(2,1,2)
levelsplot(sys,'z',[-10,10])
