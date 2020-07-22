% Imaging with calibrator transfer

for t = 1:T
    D_approx{t} = ones(4,S2,na); % Identity Jones matrices
    D_approx{t}(2,S2,:) = 0; %(0.01 + 2*A); %0;
    D_approx{t}(3,S2,:) = 0; %(0.01 + 2*A); % 0;
end

% [~,D_true] = computeD_reg(U_th, T, [], [], T); % true DDEs

parfor t = 1:T
    [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D_approx{t}, D_approx{t}, [v(:,t), u(:,t)], K, S, J, W{t}); 
%     [G1{t},G2{t},G3{t},G4{t}] = createGnufft_T2_parallel_T(D_true{t}, D_true{t}, [v(:,t), u(:,t)], K, S_true, J, W{t}); 
end

G0_{1} = cell2mat(G1_');
G0_{2} = cell2mat(G2_');
G0_{3} = cell2mat(G3_');
G0_{4} = cell2mat(G4_');

% G0{1} = cell2mat(G1'); clear G1;
% G0{2} = cell2mat(G2'); clear G2;
% G0{3} = cell2mat(G3'); clear G3;
% G0{4} = cell2mat(G4'); clear G4;

clear G1_; 
clear G2_; 
clear G3_;
clear G4_; 


for i = 1:4
    A_die{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
    At_die{i} = @(x) (so_fft2_adj((G0_{i})'*x, N, K, sp_scale));
end


Ax = 2*op_norm_stokes_RIME(A_die,At_die, size(x_th.x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
param_dual.dual_vv = 0;
param_pd.nu3 = Ax;

%%
param = param_pd;
param.evl = Ax;
for i =1:3
    epsilon{i} = zeros(size(x_th.x{i}));
    if param_pd.dual_fb
        param_l1.weights{i} = ones(size(param_l1.Psit(epsilon{i})));
        eta_{i,1} = [];
    else
        for k = 1:param_pd.Ps
            param_l1.weights{i,k} = ones(size(param_l1.Psit{k}(epsilon{i})));
            eta_{i,k} = param_l1.eta_o(i)*param_l1.weights{i,k}*1.9/Ax;
            x_th.psit_xo{i,k} = param_l1.Psit{k}(epsilon{i});
        end
    end
 end

%%
[x_approx, snr] = update_stokes(y, Jeps, tol_x, nIter, A_die, At_die,L, Lt, M, param_l1.eta_o, param_l1,param_pd,x_th, param_dual, Ax, epsilon, eta_);

for i = 1:stokes_P
    x_approx_die{i} = x_approx{i}/max(x_approx{i}(:));
    scale(i) = sum(x_approx_die{i}(:).*x_th.x{i}(:))/sum(x_approx_die{i}(:).^2);
    snr_approx_die(i) = SNR(scale(i)*x_approx_die{i}, x_th.x{i});
    snr_approx_die_unscaled(i) = snr{i}(end);
end