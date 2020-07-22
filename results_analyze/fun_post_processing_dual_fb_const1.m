function [x_final,snr] = fun_post_processing_dual_fb(A_center,A,A_off,seed,im_ch,L_ch)
%clc; clear all; close all;

addpath data
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath results
addpath ../../Tools
addpath dde_tmp_reg


T = 200;
if im_ch == 1
im_choice = 'cyg_a';
elseif im_ch == 2
im_choice = 'hydra'; %'cyg_a'; % 'avery'; % W28_256 M31
else
im_choice = 'hydra5';
end

cov_type = 'vlaa';
S = 5;
S_true = S;
data_type = 3; %1:cal_transfer; 2: cal_transfer+DDEs for off diagonal term; 3: no cal_transfer



%A_center = 5e-2; %1e-2;
%A = 5e-2; %1e-2;
%A_off = 5e-4;
test_number = 511;
off_diag = 0;
param_pd.dual_fb = 1; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox


    rng(seed);
% Loading data and initialization variables

name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_off_diag',num2str(off_diag),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed)];
n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name_ = [name, n_tmp];

name2 = name_;

D_th = [];
F =[];
J =[];
K =[];
N=[];
Omega=[];
P=[];
U_th=[];
V=[];
W=[];
Y=[];
sp_scale=[];
u=[];
v=[];
x_th=[];
y=[];
na = [];

 load(['synth_data_reg_eusipco_', name,'.mat'],'D_th','F','J','K','N','Omega','P','U_th','V','W','Y','sp_scale','u','v','x_th','y','na');
 load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die');

  
n_tmp_ = ['_no_offdiag_new_with_const_method_',num2str(param_pd.dual_fb)];
name2_ = [name, n_tmp_];

name2 = [name2_,n_tmp];

load(['joint_imaging_dde_reg_mod_pdfb_', name2_,'.mat'],'x0');
% load(['joint_imaging_dde_reg_mod_pdfb_', name2_,'.mat']);
 
% load(['results/img_dde_reg_change_',name,'_new_with_const_method_1_iter_15.mat'])
load(['img_dde_const_change_',name,'_no_offdiag_new_with_const_method_1.mat'])
%load(['joint_imaging_dde_reg_mod_pdfb_', name2_,'.mat'],'x0');

util_create_pool_bis(4, []);


param_pd.dual_fb = 1;
if param_pd.dual_fb == 0
name2 = name2_;
end


 S2 = S*S;


% --------------------------------------------------------------------------
% Parameters for primal-dual to compute the proximity operator
% --------------------------------------------------------------------------
param_pd.rel_obj = 1e-4; % stopping criterion
param_pd.max_iter = 100; % max number of iterations
param_pd.nnls_init=0;
[param_pd.Ny, param_pd.Nx] = size(x_th.x{1});
% param_pd.Nx = nnx;
% param_pd.Ny = nny;
param_pd.method = 3;
param_pd.nu1 = 1; % bound on the norm of the operator Psi
% param.nu3 = Ax %evl_precond; % bound on the norm of the operator A*G
param_pd.nu2 = 1; %for ephigraphical projection
param_pd.nu3 = 1;
param_pd.pol_thresh = 1e-2*param_pd.Ny*param_pd.Nx; % Threshold value for num of pixels not satisfying pol. const.
param_pd.pol_tol = 0; % Tolerance on the polarization constraint: \sqrt(Q^2+U^2+V^2) - I <= tol
param_pd.verbose = 4; % print log or not
param_pd.real = 1;

% -------------------------------------------------------------------------
% regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;

if param_pd.dual_fb
param_algo.opt_dict ='SARA';
[param_algo.Psi, param_algo.Psit] = SARA_sparse_operator(randn(size(x_th.x{1})), param_algo.nlevel,param_algo.opt_dict);
else
param_algo.opt_dict = {'db1', 'db2', 'db3', 'db4', 'db5', 'db6', 'db7', 'db8', 'self'};
[param_algo.Psi, param_algo.Psit] = op_p_sp_wlt_basis(param_algo.opt_dict, param_algo.nlevel, param_pd.Ny, param_pd.Nx);
end
param_pd.Ps = length(param_algo.Psit);

% -------------------------------------------------------------------------

yy = y(:);
s_col = size(yy,1)/4; 
y = reshape(yy,[s_col,4]);

param_pd.M = size(y,1); % Number of measurements for each correlator
M = param_pd.M;

if L_ch == 1
% For linear feeds
 L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
% 
% Lt = 0.5*(conj(L))'; %Adjoint conversion matrix
else
% For circular feeds
 L = [1,0,0,1;0,1,1,0;0,1i,-1i,0];
end

Lt = 0.5*L'; %0.5*(conj(L))'; %Adjoint conversion matrix

%% Joint imaging and calibration (DDEs)

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));


  
  
 %% To run dual forward backward for final imaging without polarization constraint

 %%
 [~,D1] = computeD_reg(U1, T, [], [], T); % estimated DDEs
 [~,D2] = computeD_reg(U2, T, [], [], T); % estimated DDEs
 
 parfor t = 1:T
     [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v(:,t), u(:,t)], K, S, J, W{t});
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
alpha_weights = 0.7; %* ones(1,9);
param_l1.reweight_alpha_ff = alpha_weights;
param_l1.reweight_abs_of_max = inf;
param_l1.total_reweights = 10;

util_create_pool(10);
param_dual.dual_vv = 0;

Jeps = 1500;
tol_x = 1e-4;
nIter = 1;

if strcmp(im_choice,'hydra')
    if param_pd.dual_fb
        param_l1.eta_o = [70, 20, 20]; %[43, 3, 100]; %[950, 800, 500]; %[500,400,400];
    else
        param_l1.eta_o =  [70, 20, 20]; %[550, 250, 450]; %[30,20,20]; %[500, 400, 400];
    end
else
    if param_pd.dual_fb
        param_l1.eta_o = [500, 400, 400]; %[950, 800, 500]; %[500,400,400];
    else
        param_l1.eta_o =  [30,20,20]; %[500, 400, 400];
    end
end

%param_l1.eta_o = eta_o;

eta_o = param_l1.eta_o;
for i =1:3
   x{i} = x0{i} + epsilon{i};
epsilon_{i} = zeros(size(x{i}));
if param_pd.dual_fb
   param_l1.weights{i} = ones(size(param_l1.Psit(epsilon{i})));
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

 %%
save(['results/joint_imaging_dde_reg_', name2,'_post_processing.mat'], 'x_final','eta_o');


snr_x1 = SNR(x_final{1}, x_th.x{1});
snr_x2 = SNR(x_final{2}, x_th.x{2});
snr_x3 = SNR(x_final{3}, x_th.x{3});

fprintf('%2s\t%2s\t%2s  \n', 'Final SNR_I', 'SNR_Q', 'SNR_U')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e\t%.2e\t%.2e \n',snr_x1, snr_x2, snr_x3);
fprintf('================================================================================\n');
end

