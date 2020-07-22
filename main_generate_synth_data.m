%-------------------------------------------------------------------------%
% Main file to generate the synthetic datasets for simulations.
%-------------------------------------------------------------------------%


clc; clear all; close all;
format compact;

addpath utils 
addpath utils/nufft 
addpath utils/alg
addpath lib
addpath (genpath('synth'))
addpath ../../Tools
addpath data
% addpath ../../RICK_CYG_A_Stokes_Images
% addpath ../../RICK_HYDRA_Stokes_Images


seed = 2;
rng(seed)

% Image parameters
im_choice =  'cyg_a_fits'; % Choose cyg_a or hydra5
N = 128*[1,1];
K = 2*N;
stokes_P = 3; % Number of Stokes parameters to be estimated
test_number = 1;
param_im.P = stokes_P;
param_pd.dual_fb = 0; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox


% Coverage parameters
cov_type = 'vlaa'; %choices: 'vlaa', 'meerkat'
T = 2; %200; %Number of time integrations
na = 27; % Number of antennas 
M = T*na*(na-1)/2;
n_pairs = na*(na-1)/2;


% DDE parameters
S = 5; % size of the DDE spatial support (S^2 elements in total)
J = 5; % size of the gridding kernels (J^2 elements in total)
P = 3; % size of the DDE temporal support
F = T; % size of the temporal Fourier space
dl = 1.2;  % pixel size (typically 1.5)
input_snr = 150; % noise level
A_center = 5e-2; % DDEs amplitude for diagonal terms
A = 5e-2; % DIEs amplitude for off-diagonal terms
A_off = 5e-4; % DDEs amplitude for off-diagonal terms
L_ch = 1;
off_diag = 0;

name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed),'.mat'];


% Generate data
[Y_, V_, W, Omega_, u_ab, v_ab, y, x_th, sp_scale, U_, D_, N, K, na, sig_noise] = generate_synth_data(S, J, P, T, F, dl, cov_type, A, A_center,im_choice,stokes_P,param_im,off_diag, A_off,L_ch);
x_th.x = x_th.im_true;
mkdir('data')


%%
% Save data for the calibration algorithm without temporal
% regularization
u = cell(T,1);
v = cell(T,1);
V = cell(T,1);
D_th = cell(4, T);
Omega = cell(T,1);
% y = y(:);
Y = cell(4, T);
for i = 1:4
for t = 1:T
     Y{i,t} = Y_(:, :, t,i);
    Omega{t} = Omega_(:, :, :, t);
    u{t} = u_ab(:, t);
    v{t} = v_ab(:, t);
    V{t} = V_(:, :, :, t);
    D_th{i,t} = D_{i}(:, :, t);
    U_th = U_;
end
end
save(['data/synth_data_reg_cell_', name], 'x_th', 'y', 'Y', 'Omega', 'V', 'sp_scale', 'D_th', 'U_th', 'u', 'v', 'T', 'na', 'N', 'K', 'S', 'J', 'P', 'F', 'W', 'A', 'A_center','A_off','sig_noise');

% % Rename the data for calibration with temporal regularization
u = u_ab;
v = v_ab;
V = V_;
D_th = D_;
U_th = U_;
Omega = Omega_;
Y = Y_;
y = y(:);
save(['data/synth_data_reg_', name], 'x_th', 'y', 'Y', 'Omega', 'V', 'sp_scale', 'D_th', 'U_th', 'u', 'v', 'T', 'na', 'N', 'K', 'S', 'J', 'P', 'F', 'W', 'A', 'A_center','A_off','sig_noise');


