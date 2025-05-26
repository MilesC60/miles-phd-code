%% yeah lmao
% ignore all refs to S_1, I took it out bs we don't need it
% This function solves for triplet-pair dynamics using SLE_fast.
% DM stands for density matrix, low-J and high-J are called conformers sometimes
function [rhot] = fluctuating_J_populations(X,model,ks,ts)
% X is the vector of all 15 spin Ham parameters, in this order:
if numel(X) ~= 15
    error('Wrong no. spin Ham parameters! Should be 12: [D1,E1,D2,E2,Jlow,Jhigh,eul1,eul2,bvec]')
end
% D/E are ZFS parameters for triplets 1 and 2, Jlow and Jhigh are the J-values
% between which we fluctuate, eul1, eul2, bvec are 1x3 arrays (Euler angles for
% each triplet, and field direction)
% bvecs is a 3xNB list of field vectors (Bx,By,Bz)

% model is one of three models: "pulse" means that the system starts in one state (e.g 1TT)
% with no inflowing population; "cw" means the initial state is zero, and population
% is inflowing; "steady" is like "cw" but just computes the steady state

% ks is the list of rates in the following order:
if numel(ks) ~= 2
    error('Wrong no. rates! Should be 2: [k_tt,k_flap]')
end
% TT decay rate; J-fluctuation frequency;

% pull out spin Ham parameters from X
D1=X(1);E1=X(2);D2=X(3);E2=X(4);Jlow = X(5);Jhigh = X(6);
eul1 = X(7:9);eul2 = X(10:12);bvec = X(13:15);
% read rates from the function input 
k_tt = ks(1);k_flap = ks(2);

% Define spin Ham
sys.S = [1 1]; sys.J = 1;% just to get the form of isotropic exchange Hamiltonian
[Hexc,mux,muy,muz] = ham(sys);% Get Hexc and parts of Zeeman Hamiltonian
sys.J=0;% so we can get the ZFS Hamiltonian
sys.D = [D1,E1;D2,E2]; sys.DFrame = [eul1;eul2];
[Hzfs,~,~,~] = ham(sys);% get ZFS ham
Hzee =  -(mux*bvec(1) + muy*bvec(2) + muz*bvec(3));
H_both = Hzee + Hzfs;% this is everything but exchange

h_high = Jhigh*Hexc + H_both;% high-J Hamiltonian
h_low = Jlow*Hexc + H_both;% low-J Ham
[v_high,~] = eig(h_high);[v_low,~] = eig(h_low);% get the high-J and low-J eigenstates
% to later convert the DMs into the correct basis

% The total Hamiltonian is 19x19; 1 S_1 + 9 TT at high J + 9 TT at low J
% diagonally stack 0 then H_high then H_low, i.e S_1 has zero energy (doesn't matter)
H = blkdiag(h_high,h_low);

ss = cgmatrix(1,1);tt1 = ss(9,:)';% tt1 is pure 1TT, the S^2 eigenstate
tt1h = [tt1;zeros(9,1)];tt1l = [zeros(9,1);tt1];% get 1TT in both conformers
rho0 = (tt1h*tt1h' + tt1l*tt1l')/2;

% Define the jump operators
J_up = kron([0 1;0 0],eye(9));% from low-J to high-J conformer
J_down = J_up';% J_down is J_up in reverse

j_list = {J_up,J_down};% make the lists of operators and rates
jk_list = {k_flap,k_flap};% could change fourth rate for detailed balance

% Define decay operators
K_tt = eye(18);% everything else

k_list = {K_tt};% make the lists of operators and rates
kk_list = {k_tt};

rhot = [];% empty array for TT pops
rhot_big = SLE_code(H,j_list,jk_list,k_list,kk_list,rho0,ts,model);
% rhot_big returns a 18-by-(18xlength(ts)) array, DMs stacked next to each other
for m = 1:length(ts)
    %rhot_rho = reshape(rhot_big(:,m),18,18);
    rho = reshape(rhot_big(:,m),18,18);
    %rho = rhot_rho(:,18*(m-1)+1:18*m);% select one 19-by-19 DM

    rho_high = rho(1:9,1:9);% take the high-J 9x9 DM
    p_high = v_high'*rho_high*v_high;% convert it to the H eigenbasis

    rho_low = rho(10:18,10:18);% repeat for the low-J 9x9 DM
    p_low = v_high'*rho_low*v_high;
    % add average of high J and low J DMs to the output array
    rhot = [rhot,p_high+p_low];
end

end
























