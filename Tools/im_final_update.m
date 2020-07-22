% Script for final imaging update step with polarization constraint


 [~,D1] = computeD_reg(U1, T, [], [], T); % estimated DDEs
 [~,D2] = computeD_reg(U2, T, [], [], T); % estimated DDEs
 
 parfor t = 1:T
     [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v(:,t), u(:,t)], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
 end
 G0_{1} = cell2mat(G1_');
 G0_{2} = cell2mat(G2_');
 G0_{3} = cell2mat(G3_');
 G0_{4} = cell2mat(G4_');
 
 clear G1_
 clear G2_
 clear G3_
 clear G4_
 
 for i = 1:4
     A_est{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
     At_est{i} = @(x) (so_fft2_adj((G0_{i})'*x, N, K, sp_scale));
 end

  
Ax = 2*op_norm_stokes_RIME(A_est,At_est, size(x_th.x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
param.evl = Ax;

param_pd.dual_fb = 0;
param_algo.opt_dict = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
[param_algo.Psi, param_algo.Psit] = op_p_sp_wlt_basis(param_algo.opt_dict, param_algo.nlevel, param_pd.Ny, param_pd.Nx);
param_pd.Ps = length(param_algo.Psit);


%%
param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.approx = 0; % If box bounds on the image to be considered
if param_l1.approx == 0
    for i=1:3
        param_l1.xA{i} = zeros(size(x_th.x{i}));
    end 
end
param_l1.reweight_rel_obj = 5e-5; % criterion for performing reweighting
param_l1.reweight_alpha = 1; %ones(1,9) ;
alpha_weights = 0.7; % * ones(1,9);
param_l1.reweight_alpha_ff = alpha_weights;
param_l1.reweight_abs_of_max = inf;
param_l1.total_reweights = 10;

util_create_pool(param_pd.Ps);
param_dual.dual_vv = 0;

Jeps = 1500;
tol_x = 1e-4;
nIter = 1;
 if strcmp(im_choice,'hydra5')
param_l1.eta_o = [70, 20, 20];
else
            param_l1.eta_o = [30, 20, 20];
end
                                   
                                   
eta_o = param_l1.eta_o;
for i =1:3
    x{i} = x0{i} + epsilon{i};
    epsilon_{i} = zeros(size(x_th.x{i}));
    if param_pd.dual_fb
        param_l1.weights{i} = ones(size(param_l1.Psit(epsilon_{i})));
        eta_{i,1} = [];
    else
        for k = 1:param_pd.Ps
            param_l1.weights{i,k} = ones(size(param_l1.Psit{k}(epsilon_{i})));
            eta_{i,k} = param_l1.eta_o(i)*param_l1.weights{i,k}*1.9/Ax;
            x_th.psit_xo{i,k} = param_l1.Psit{k}(epsilon_{i});

        end
    end
end

[x_final, snr] = update_stokes(y, Jeps, tol_x, nIter, A_est, At_est,L, Lt, M, param_l1.eta_o, param_l1,param_pd,x_th, param_dual, Ax, x, eta_);



save(['results/joint_imaging_dde_reg_', name2,'_final.mat'], 'x_final','eta_o');


snr_x1 = SNR(x_final{1}, x_th.x{1});
snr_x2 = SNR(x_final{2}, x_th.x{2});
snr_x3 = SNR(x_final{3}, x_th.x{3});

fprintf('%2s\t%2s\t%2s  \n', 'Final SNR_I', 'SNR_Q', 'SNR_U')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e\t%.2e\t%.2e \n',snr_x1, snr_x2, snr_x3);
fprintf('================================================================================\n');

