%-------------------------------------------------------------------------%
% Main file to obtain the reference imaging results (with the calibrator
% transfer).
%-------------------------------------------------------------------------%
%% 

addpath data
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath Tools


im_choice = 'cyg_a_fits'; % Choose cyg_a or hydra5
cov_type = 'vlaa'; 
stokes_P = 3; % Number of Stokes parameters to be estimated
test_number = 1;
A_center = 5e-2; 
A = 5e-2; 
A_off = 5e-4;
T  = 2;
off_diag = 0;
L_ch = 1; %1: Linear feeds
seed = 2;
rng(seed);
m_ch = 0; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox
S = 5; 
S_true = S;

% Data loading
name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed)];

load(['synth_data_reg_', name,'.mat'],'D_th','F','J','K','N','Omega','P','U_th','V','W','Y','sp_scale','u','v','x_th','y','na');


%% Define parameters for imaging 

parameters_script; % Define parameters

% Assuming calibrator transfer is performed; zero off-diagonal terms
S = 1;
S2 = S.^2;

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));


%%
% Imaging with calibrator transfer

imaging_cal_transfer;

% x_approx: Approximated images
% snr: SNR of approximated images

disp(['SNR(x_approx_die_unscaled) = ', num2str(snr_approx_die_unscaled(1)),'  ',num2str(snr_approx_die_unscaled(2)),'  ',num2str(snr_approx_die_unscaled(3))]);
disp(['SNR(x_approx_die) = ', num2str(snr_approx_die(1)),'  ',num2str(snr_approx_die(2)),'  ',num2str(snr_approx_die(3))]);

mkdir('results')
save(['results/imaging_norm_die_', name2, '.mat'], 'eta_o', ...
       'x_approx',  'snr_approx_die_unscaled', 'x_approx_die', 'snr_approx_die','L');
   
