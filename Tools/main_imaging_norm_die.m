%-------------------------------------------------------------------------%
% Main file to obtain the reference imaging results (with the true DDEs, 
% and normalized DIEs).
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

clc; clear all; close all;

addpath data
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath ../Tools

T  = 200;
im_choice = 'cyg_a'; % 'avery'; %'cyg_a'; % 'CYGC_alex-1024'; %'cyg_a'; %'ps_256'; %'M31'; %'rand_im'; %'3c286'; %'cyg_a'; %'avery' %'rand_im'; %'M31'; %rand_im: to generate point sources image
cov_type = 'vlaa' %'random' %'vlaa';
test_aux = '';
stokes_P = 3; % Number of Stokes parameters to be estimated
data_type = 3; %1:cal_transfer; 2: cal_transfer+DDEs for off diagonal term; 3: no cal_transfer
A_center = 5e-2; %5e-2;
A = 5e-2;
test_number = 1;
off_diag = 0;
rng(2);
%  util_create_pool_bis(4, []);

S = 5;
S_true = S;

% Data loading
% name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_ data_type',num2str(data_type)]; %,'_off_diag',num2str(off_diag)];
% name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_off_diag',num2str(off_diag),'_A_center=',num2str(A_center),'_A=',num2str(A)];
name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_off_diag',num2str(off_diag),'_A_center=',num2str(A_center),'_A=',num2str(A),'_data_type=',num2str(data_type)];


load(['synth_data_reg_eusipco_', name, '.mat']) % for temporal regularization
% load(['synth_data_reg_cell_',name2,'.mat']) % without temporal regularization

[nny,nnx] = size(x_th.im_true{1});
 Ny = nny;
 Nx = nnx;
 
%% Imaging with the true DDEs
% --------------------------------------------------------------------------
% Parameters for primal-dual to compute the proximity operator
% --------------------------------------------------------------------------
param_pd.rel_obj = 1e-4; % stopping criterion
param_pd.max_iter = 200; % max number of iterations
param_pd.nnls_init=0;
[param_pd.Ny, param_pd.Nx] = size(x_th.x{1});
% param_pd.Nx = nnx;
% param_pd.Ny = nny;
param_pd.method = 3;
param_pd.nu1 = 1; % bound on the norm of the operator Psi
% param.nu3 = Ax %evl_precond; % bound on the norm of the operator A*G
param_pd.nu2 = 1; %for ephigraphical projection
param_pd.pol_thresh = 1e-3*param_pd.Ny*param_pd.Nx; % Threshold value for num of pixels not satisfying pol. const.
param_pd.pol_tol = 0; % Tolerance on the polarization constraint: \sqrt(Q^2+U^2+V^2) - I <= tol
param_pd.verbose = 4; % print log or not
param_pd.dual_fb = 0; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox
param_pd.real = 1;
param_pd.nnls_init = 0;


n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp];


% -------------------------------------------------------------------------
% regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;

if param_pd.dual_fb
param_algo.opt_dict = 'SARA'; % 'DB'; % 'SARA'; %'Dirac'; %'DB8'; % 'SARA';
[param_algo.Psi, param_algo.Psit] = SARA_sparse_operator(randn(size(x_th.x{1})), param_algo.nlevel,param_algo.opt_dict);
else
param_algo.opt_dict = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
[param_algo.Psi, param_algo.Psit] = op_p_sp_wlt_basis(param_algo.opt_dict, param_algo.nlevel, param_pd.Ny, param_pd.Nx);

for i =1:3
    for k =1:length(param_algo.Psit)
        x_th.psit_xo{i,k} = zeros(size(param_algo.Psit{k}(x_th.x{i})));
    end
end
end
param_pd.Ps = length(param_algo.Psit);

param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
% param_l1.weights = 1;
param_l1.eta_o = [1e-4, 1e-3, 1e-3];
param_l1.approx = 0; % If box bounds on the image to be considered
if param_l1.approx == 0
    for i=1:3
        param_l1.xA{i} = zeros(size(x_th.x{i}));
    end 
end
% param_l1.reweight_steps = [500:step:last_reweight  inf];
param_l1.reweight_rel_obj = 5e-5; % criterion for performing reweighting
% param_l1.reweight_min_steps_rel_obj = step;
% param_l1.reweight_max_reweight_itr = max(500,param_l1.max_iter - 500);
param_l1.reweight_alpha = 1; %ones(1,9) ;
alpha_weights = 0.7; % * ones(1,9);
param_l1.reweight_alpha_ff = alpha_weights;
param_l1.reweight_abs_of_max = inf;
param_l1.total_reweights = 10;

util_create_pool(param_pd.Ps);

% -------------------------------------------------------------------------

yy = y(:);
s_col = size(yy,1)/4; 
y = reshape(yy,[s_col,4]);

param_pd.M = size(y,1); % Number of measurements for each correlator
M = param_pd.M;

L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix

Lt = 0.5*(conj(L))'; %Adjoint conversion matrix


% -------------------------------------------------------------------------


 %% Results with approximate DDEs
 
 % Assuming calibration transfer is performed; zero off-diagonal terms
S = 1;
S2 = S.^2;
for t = 1:T
    D_approx{t} = ones(4,S2,na);
     D_approx{t}(2,S2,:) = 0; %(0.01 + 2*A); %0;
     D_approx{t}(3,S2,:) = 0; %(0.01 + 2*A); % 0;
    
end


[~,D_true] = computeD_reg(U_th, T, [], [], T); % true DDEs

% for i = 1:4
%     for t=1:T
%    D_true{i,t} = D_th{i}(:, :, t);
%     end
% end


parfor t = 1:T
          [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D_approx{t}, D_approx{t}, [v(:,t), u(:,t)], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
          [G1{t},G2{t},G3{t},G4{t}] = createGnufft_T2_parallel_T(D_true{t}, D_true{t}, [v(:,t), u(:,t)], K, S_true, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
end

G0_{1} = cell2mat(G1_');
G0_{2} = cell2mat(G2_');
G0_{3} = cell2mat(G3_');
G0_{4} = cell2mat(G4_');

G0{1} = cell2mat(G1');
G0{2} = cell2mat(G2');
G0{3} = cell2mat(G3');
G0{4} = cell2mat(G4');

clear G1_; clear G1
clear G2_; clear G2
clear G3_; clear G3
clear G4_; clear G4


for i = 1:4
    A_die{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
    At_die{i} = @(x) real(so_fft2_adj((G0_{i})'*x, N, K, sp_scale));
    
    A_true{i} = @(x) G0{i}*so_fft2(x, K, sp_scale);
    At_true{i} = @(x) real(so_fft2_adj((G0{i})'*x, N, K, sp_scale));
end

Ax = 2*op_norm_stokes_RIME(A_true,At_true, size(x_th.x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));
param_dual.dual_vv = 0;
param_pd.nu3 = Ax;

%%
param = param_pd;
param.evl = Ax;
% pdfb_stokes_imaging_true_ddes;

%%
 Jeps = 1000;
tol_x = 1e-4;
nIter = 1;
param_l1.eta_o =  [30,20,20]; %[500, 400, 400];  
eta_o = param_l1.eta_o;
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
        %     eta{i} = param_l1.eta_o(i)*weights{i}*1.9/Ax;
    end
end

% Imaging with true DDEs - for reference 
[xq, snr] = update_stokes(y, Jeps, tol_x, nIter, A_true, At_true,L, Lt, M, param_l1.eta_o, param_l1,param_pd,x_th, param_dual, Ax, epsilon, eta_);

for i = 1:stokes_P
x_true_dde{i} = xq{i}/max(xq{i}(:));
    scale(i) = sum(x_true_dde{i}(:).*x_th.x{i}(:))/sum(x_true_dde{i}(:).^2);
    snr_true_dde(i) = SNR(scale(i)*x_true_dde{i}, x_th.x{i});
    snr_true_dde_unscaled(i) = snr{i}(end);
end

disp(['SNR(x_true_dde_unscaled) = ', num2str(snr_true_dde_unscaled(1)),'  ',num2str(snr_true_dde_unscaled(2)),'  ',num2str(snr_true_dde_unscaled(3))]);
disp(['SNR(x_true_dde) = ', num2str(snr_true_dde(1)),'  ',num2str(snr_true_dde(2)),'  ',num2str(snr_true_dde(3))]);

%%
% Imaging with normalized DIEs
Ax = 2*op_norm_stokes_RIME(A_die,At_die, size(x_th.x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
param.evl = Ax;
%%
[x_approx, snr] = update_stokes(y, Jeps, tol_x, nIter, A_die, At_die,L, Lt, M, param_l1.eta_o, param_l1,param_pd,x_th, param_dual, Ax, epsilon, eta_);

for i = 1:stokes_P
x_approx_die{i} = x_approx{i}/max(x_approx{i}(:));
scale(i) = sum(x_approx_die{i}(:).*x_th.x{i}(:))/sum(x_approx_die{i}(:).^2);
snr_approx_die(i) = SNR(scale(i)*x_approx_die{i}, x_th.x{i});
snr_approx_die_unscaled(i) = snr{i}(end);
end

disp(['SNR(x_approx_die_unscaled) = ', num2str(snr_approx_die_unscaled(1)),'  ',num2str(snr_approx_die_unscaled(2)),'  ',num2str(snr_approx_die_unscaled(3))]);
disp(['SNR(x_approx_die) = ', num2str(snr_approx_die(1)),'  ',num2str(snr_approx_die(2)),'  ',num2str(snr_approx_die(3))]);

mkdir('results')
save(['results/imaging_norm_die_', name2, '.mat'], 'xq', 'snr_true_dde_unscaled', 'eta_o', ...
      'x_true_dde', 'snr_true_dde', 'x_approx',  'snr_approx_die_unscaled', 'x_approx_die', 'snr_approx_die');
  
  
