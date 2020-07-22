% function [x_final,snr] = fun_main_imaging_dde_reg_new1(A_center,A,A_off,seed,im_ch,L_ch)

%-------------------------------------------------------------------------%
% Main file to run joint calibration and imaging algorithm.
%-------------------------------------------------------------------------%
%% 

addpath data
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath results
addpath Tools
addpath dde_tmp_reg


im_choice = 'cyg_a_fits'; % Choose cyg_a or hydra5
T = 200;
cov_type = 'vlaa'; % vlaa, meerkat
S = 5; % DDE support size in spatial Fourier domain
S_true = S;
seed = 2;
rng(seed);

% Parameters to constrain amplitudes of DIE/DDE Fourier kernels
A_center = 5e-2; % DDEs amplitude for diagonal terms
A = 5e-2; % DIEs amplitude for off-diagonal terms
A_off = 5e-4; % DDEs amplitude for off-diagonal terms
L_ch = 1;
off_diag = 0;
m_ch = 0; % 0: run primal dual to compute prox; 1: run dual forward backward to compute prox
param_pd.dual_fb = m_ch;

% Loading data and initialization variables
name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed)];

n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name_ = [name, n_tmp];

load(['synth_data_reg_', name,'.mat'],'D_th','F','J','K','N','Omega','P','U_th','V','W','Y','sp_scale','u','v','x_th','y','na');

load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','x_approx_die','snr_approx_die')
util_create_pool_bis(12, []);

S2 = S*S;

n_tmp_ = ['_new_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];


%% Define parameters

parameters_joint_cal_im;

%% Run joint calibration and imaging algorithm

[U1, U2, epsilon, objective, time_dde, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, ...
    U2, sp_scale, V, param_im, param_dde, param_algo, W, name2, mu,param_pd,L,Lt,eta__,U_th);

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
 

