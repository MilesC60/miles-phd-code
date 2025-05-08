%% 
clear
% use EasySpin to set up the matrices
sys.S = [1 1];
sys.D = [1414,14;1414,14];
sys.J = 0;
[h0,mux,muy,muz] = ham(sys);
bvec = [0,0,1];bvec = bvec/norm(bvec);
bvals = -200:1:200;nb = numel(bvals);
pops = zeros(1,nb);
k_d = 1e2;

ss = cgmatrix(1,1); tt1 = ss(9,:)';
vec_tt1 = reshape(tt1*tt1',[],1);
rho0 = tt1*tt1';
vec_rho0 = reshape(rho0,[],1);
for m = 1:nb
    hzee = bvals(m)*mux;
    h = h0 + hzee;
    H = -1j*2*pi*(kron(eye(9),h)-kron(transpose(h),eye(9)));
    K = k_d*eye(81);
    A = H - K;

    vec_ss = -A\vec_rho0;
    pops(m) = real(vec_ss'*vec_tt1);
end

figure('Position',[89 500 746 318])
plot(bvals,pops,LineWidth=2)
xlabel('Field Strength (mT)')
ylabel('MPL Intensity (a.u.)')
