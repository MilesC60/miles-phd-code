% This function solves the Stochastic Liouville Equation in one of three models:
% "pulse" means that the system starts in one state, with no inflowing
% population; "cw" means the initial state is zero, and population is
% inflowing; "steady" is like "cw" but just computes the steady state
function [vec_rho_t] = SLE_fast(H,j_list,jk_list,k_list,kk_list,rho0,ts,model)
% h is the system Hamiltonian, which has to be time-independent sorry
% j_list is a cell list {j1,j2,...} of jump operators; e.g something like
% j1 = |a><b| drives incoherent transitions from b to a.
% jk_list is a cell list {k1,k2,...} that tells the jump operators in
% j_list how fast to go. The order of jk_list has to match j_list, i.e the
% rate jk_list{1} tells j_list{1} how fast to go
% k_list is a cell list of decay operators, e.g |a><a| causes the a
% population to decay
% kk_list is a cell list of decay rates
% rho0 is the inital state ("pulse") or inflowing state ("cw" or "steady")
% ts is a list of timesteps, must be a singleton for "steady"
% model is the model, "steady" or "cw" or "pulse"

% Catch some common errors
if numel(j_list) ~= numel(jk_list)
    error('Wrong no. jump rates and operators!')
end
if numel(k_list) ~= numel(kk_list)
    error('Wrong no. decay rates and operators!')
end
if (model == "steady") & (length(ts) ~= 1)
    error('Please set ts to a list with length 1 (or choose non steady-state!)')
end

N = size(rho0,1);I = eye(N);% N is the Hilbert space dimension
vec_rho_0 = reshape(rho0,[],1);% vectorise (stack columns) rho0
vec_rho_t = zeros(length(vec_rho_0), length(ts));% array to be filled

hbar = 1/2/pi;% so time in microseconds, energy in MHz
L = -(1i/hbar)*(kron(I,H)-kron(H.',I));% initialise the Liouville operator,
% we will add terms to this as we go

% add all decay terms
for m = 1:length(kk_list)
    K = k_list{m};
    gamma = kk_list{m};
    L = L - gamma*0.5*(kron(I,K)+kron(transpose(K),I));
end
% add all jump operator terms
for m = 1:length(jk_list)
    J = j_list{m};
    gamma = jk_list{m};
    L = L + gamma*(kron(conj(J), J) - (1/2)*(kron((J'*J).',I) + kron(I, J'*J)));
end

% Computing stuff with the Liouville operator is quicker with its
% eigendecomposition, but since it's not Hermitian or unitary you have to
% use left/right eigenvectors. All the stuff using R_mat, D, L_mat is
% copied from Francesco's OQS tutorial 10.1103/PRXQuantum.5.020202 

% In "cw" model I assume there's a steady state, otherwise you get
% exponential growth (and I have to work harder to solve it)


if model == "steady"
    vec_rho_t(:,1) = -L\vec_rho_0;% just -(L^-1)b, c.f. my thesis
else
    [R_mat, D, ~] = eig(L);
    if model == "cw"
        if min(abs(diag(D))) < 1e-10% this means one state doesn't decay
            error("No steady state! Add more decay rates or remove states that don't decay")
        end
    end
    R = [R_mat vec_rho_0];% Concatenating R_mat and vec_rho_0
    R = rref(R);  
    %tbh just read Francesco's paper for this
    % for k=1:(N^2)
    %     b(k) = R(k,N^2+1);
    %     if model == "pulse"% then rho0 is the INITIAL density matrix, no incoming DM
    %         vec_rho_t(:,:) = vec_rho_t(:,:) + b(k)*R_mat(:,k)*exp(D(k,k)*ts);
    %         % this is exp(Lt)
    %     elseif model == "cw"% then rho0 is the INCOMING density matrix, no initial DM
    %         vec_rho_t(:,:) = vec_rho_t(:,:) + b(k)*R_mat(:,k)*(exp(D(k,k)*ts)-1)/D(k,k);
    %         % this is (exp(Lt)-I)L^-1
    %     else
    %         error("I don't understand "+'"'+model+'"')
    %     end
    % end
    if model == "pulse"
        bd = R(:,N^2+1).*(exp(diag(D)*ts));
        vec_rho_t = vec_rho_t + R_mat*bd;
    elseif model == "cw"
        bd = R(:,N^2+1).*((exp(diag(D)*ts)-1)./diag(D));
        vec_rho_t = vec_rho_t + R_mat*bd;
    end
end

end




















