% function [x_final,snr] = fun_main_imaging_dde_reg_new1(A_center,A,A_off,seed,im_ch,L_ch)

% clc; clear all; close all;

addpath(genpath('../../Data/jbirdi/Pol_cal/final_cirrus'))
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath results
addpath Tools
addpath dde_tmp_reg
addpath init_scripts


im_choice = 'cyg_a'; % Choose image: cyg_a or hydra5
T = 200;

im_ch = 1;

if im_ch == 1
im_choice = 'cyg_a';
elseif im_ch == 2
im_choice = 'hydra'; %'cyg_a'; % 'avery'; % W28_256 M31
else
im_choice = 'hydra5';
end

cov_type = 'vlaa'; % vlaa, meerkat
S = 5; % DDE support size in spatial Fourier domain
S_true = S;
seed = 2;

rng(seed);
test_number = 11;

% Parameters to constrain amplitudes of DDE Fourier kernels
A_center = 5e-2; % variance for DIEs of diagonal Jones terms
A = 5e-2; % variance for DIEs of off-diagonal Jones terms 
A_off = 5e-4; % variance for DDEs of off-diagonal Jones terms

L_ch = 1; % 1 for linear feeds; 2 for circular feeds


off_diag = 0;

% Parameters for algorithm computing proximity operator of the image
% regularization term
param_pd.dual_fb = 1; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox

% Loading data and initialization variables
name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_off_diag',num2str(off_diag),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed)];

n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name_ = [name, n_tmp];

load(['synth_data_reg_eusipco_', name,'.mat'],'D_th','F','J','K','N','Omega','P','U_th','V','W','Y','sp_scale','u','v','x_th','y','na');

load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die')
util_create_pool_bis(12, []);

S2 = S*S;

n_tmp_ = ['_new_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];


% -------------------------------------------------------------------------
% regularisation parameters
% -------------------------------------------------------------------------
param_algo.eta = eta_o; % Regularization parameters for the images
param_algo.Jeps = 100; % Number of FB iterations for imaging
param_algo.tol_x = 1e-5; % Stopping criterion for imaging
param_algo.tol_crit = 2e-2; % Global stopping criterion
% parameters for DDE estimation
param_dde.max_it = 20; % Global number of iterations
param_dde.JUtot = 0; %2; % Number of iterations to complete a cycle
param_dde.JU1 = 5;   % Number of FB iterations to update U1
param_dde.JU2 = 5;   % Number of FB iterations to update U2
param_dde.nu = 1000; % Regularization parameter for DDEs update
param_dde.tol_norm = 1e-5; % Stopping criterion for DDEs update
param_dde.S = S;
param_dde.A_center = A_center;
param_dde.A = A;



% --------------------------------------------------------------------------
% Parameters for primal-dual to compute the proximity operator
% --------------------------------------------------------------------------
param_pd.rel_obj = 1e-4; % stopping criterion
param_pd.max_iter = 100; % max number of iterations
param_pd.nnls_init=0;
[param_pd.Ny, param_pd.Nx] = size(x_th.x{1});
param_pd.method = 3;
param_pd.nu1 = 1; % bound on the norm of the operator Psi
% param.nu3 = Ax %evl_precond; % bound on the norm of the operator A*G
param_pd.nu2 = 1; %for ephigraphical projection
param_pd.nu3 = 1;
param_pd.pol_thresh = 1e-2*param_pd.Ny*param_pd.Nx; % Threshold value for num of pixels not satisfying pol. const.
param_pd.verbose = 4; % print log or not
param_pd.real = 1;


% -------------------------------------------------------------------------
% percentage for the threshold of the image after reconstruction with DIEs
% -------------------------------------------------------------------------
param_algo.Tau = 0.01;
param_algo.p_min = 0.2;
param_algo.p_max = 1.1;
% dwtmode('per'); % per

% -------------------------------------------------------------------------
% regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;
param_pd.dual_fb = 0;


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
%  L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
 L = [1,0,0,1;1,0,0,-1]; % Conversion matrix
else
% For circular feeds
 L = [1,0,0,1;0,1,1,0;0,1i,-1i,0];
end
Lt = 0.5*L'; %Adjoint conversion matrix


%% Joint imaging and calibration (DDEs)

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

ddes_init;
im_init;

%%
% DDE estimation (set parameters & estimate DDE + image)


B_tmp = 1;
Bt_tmp = 1;

if strcmp(im_choice,'hydra5')
   if param_pd.dual_fb
    param_l1.eta_o = [70, 20, 20];
else
param_l1.eta_o =  [70, 20, 20];
   end
else
if param_pd.dual_fb
    param_l1.eta_o = [500, 400, 400]; 
else
param_l1.eta_o =  [30,20,20];
end
end
%%
param_l1.eta_o = [500, 400, 400]; 
param_pd.dual_fb = 2;

%%

eta__ = param_l1.eta_o;

param_algo.eta = eta__;
[U1, U2, epsilon, objective, time_dde, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate_no_U(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, ...
    U2, sp_scale, V, param_im, param_dde, param_algo, W, name2,param_pd,L,Lt,eta__,U_th);



%%
for i =1:3
    x{i} = x0{i} + epsilon{i};
    x_reg{i} = x{i}/max(x{i}(:));
    scale_x(i) = sum(x_reg{i}(:).*x_th.x{i}(:))/sum(x_reg{i}(:).^2);
    snr_reg(i) = SNR(scale_x(i)*x_reg{i}, x_th.x{i});
end

mkdir('results')
save(['results/joint_imaging_dde_reg_', name2,'.mat'], 'U1', 'U2', 'epsilon', ...
    'x0', 'objective', 'param_dde', 'time_dde', 'error_dirty', 'snr_dirty', 'snr_x', ...
    'snr_reg', 'x_reg', 'snr_approx_die', 'x_approx_die','param_algo','param_pd','L');



  
 %% Final imaging step with estimated DDEs
 
 im_final_update;
 

