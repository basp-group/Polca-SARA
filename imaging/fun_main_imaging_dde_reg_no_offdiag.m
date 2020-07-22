function [x_final,snr] = fun_main_imaging_dde_reg_no_offdiag(A_center,A,A_off,seed,im_ch,L_ch)

% clc; clear all; close all;

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

cov_type = 'vlaa'; % vlaa, meerkat
test_aux = '';
S = 5;
S_true = S;
data_type = 3; %1:cal_transfer; 2: cal_transfer+DDEs for off diagonal term; 3: no cal_transfer
% seed = 2;
rng(seed);
% A_center = 5e-2; %5e-2; %5e-2;
% A = 5e-2; %5e-2;
% A_off = 5e-4;
if L_ch == 1
test_number = 411; %211;
else
test_number = 611;
end

if strcmp(im_choice,'hydra5')
test_number = 511;
end

off_diag = 0;
param_pd.dual_fb = 0; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox

% Loading data and initialization variables
name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_off_diag',num2str(off_diag),'_A_center=',num2str(A_center),'_A=',num2str(A),'_seed=',num2str(seed)];

n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name_ = [name, n_tmp];

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

% name = [cov_type, '_', im_choice, '_', test_number,'_',num2str(S)];
% load(['synth_data_reg_eusipco_', name,'.mat'])
% load(['synth_data_reg_cell_', name, '.mat'])
 load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die')
util_create_pool_bis(12, []);

S2 = S*S;

n_tmp_ = ['_no_off_diag_new_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];


% -------------------------------------------------------------------------
% regularisation parameters
% -------------------------------------------------------------------------
param_algo.eta = eta_o; %1e-3; %eta;
param_algo.Jeps = 100; %100; % 1000
param_algo.tol_x = 1e-5;
param_algo.tol_crit = 2e-2; %2e-2;
% parameters for DDE estimation
param_dde.max_it = 20; %20;
param_dde.JUtot = 2; %5; % 10;
param_dde.JU1 = 5;   % 5
param_dde.JU2 = 5;   % 5
param_dde.nu = 1000; % 1000;
mu = 0;
param_dde.tol_norm = 1e-5;
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
% percentage for the threshold of the image after reconstruction with DIEs
% -------------------------------------------------------------------------
param_algo.Tau = 0.01;
param_algo.p_min = 0.2;
param_algo.p_max = 1.1;
% dwtmode('per'); % per
Jt = [];

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


c = floor(S2/2) + 1;
p = floor(P/2) + 1;

% DDE parameters
reg = 1; % 0: no temporal regularization; 1: temporal regularization
if reg == 0 

param_dde.theta_maxR = 2*A_center*ones(S2, P, 4); %*sqrt(F); 
param_dde.theta_maxR(:,:,2) = 0;
param_dde.theta_maxR(:,:,3) = 0;
param_dde.theta_maxR(c,:,2) = 2*A*ones(1,P,1); %*sqrt(F);
param_dde.theta_maxR(c,:,3) = 2*A*ones(1,P,1); %*sqrt(F);
param_dde.theta_minR = -param_dde.theta_maxR;
param_dde.theta_maxI = param_dde.theta_maxR;
param_dde.theta_minI = -param_dde.theta_maxI;

% For DIEs
for i =1:3:4
param_dde.theta_maxR(c,p,i) = (1 + 2*A_center); %*sqrt(F);
param_dde.theta_minR(c,p,i) = (1 - 2*A_center); %*sqrt(F);
% param_dde.theta_maxI(c,p,i) = 0; %2*A*ones(S2, P);
% param_dde.theta_minI(c,p,i) = -param_dde.theta_maxI(c,p,i);
end
for i = 2:3
param_dde.theta_maxR(c,p,i) = (off_diag + 2*A); %*sqrt(F); % 1 + 2*die_var
% param_dde.theta_maxI(c,p,i) = 0; %(0.1+2*A); %2*A_center;       % 2*die_var
param_dde.theta_minR(c,p,i) = (off_diag - 2*A); %*sqrt(F); % 1 - 2*die_var
% param_dde.theta_minI(c,p,i) = 0; %-(0.1-2*A); %-2*A_center;
end

else

param_dde.theta_maxR = A_center*ones(S2, P, 4)*sqrt(F);
param_dde.theta_maxR(:,:,2) = 0; %A_off*ones(S2, P, 1)*sqrt(F);  %0;
param_dde.theta_maxR(:,:,3) = 0; %A_off*ones(S2, P, 1)*sqrt(F); %0;
param_dde.theta_maxR(:,p,:) = 0;
%param_dde.theta_maxR(c,:,2) = A*ones(1,P,1)*sqrt(F);
%param_dde.theta_maxR(c,:,3) = A*ones(1,P,1)*sqrt(F);
param_dde.theta_minR = -param_dde.theta_maxR;
param_dde.theta_maxI = param_dde.theta_maxR;
param_dde.theta_minI = -param_dde.theta_maxI;

% For DIEs
for i =1:3:4
param_dde.theta_maxR(c,p,i) = sqrt(F);
param_dde.theta_minR(c,p,i) = sqrt(F);
param_dde.theta_maxI(c,p,i) = 0; %2*A*ones(S2, P);
param_dde.theta_minI(c,p,i) = -param_dde.theta_maxI(c,p,i);
end
%for i = 2:3
%param_dde.theta_maxR(c,p,i) = (off_diag)*sqrt(F); % 1 + 2*die_var
%param_dde.theta_maxI(c,p,i) = 0; %(0.1+2*A); %2*A_center;       % 2*die_var
%param_dde.theta_minR(c,p,i) = (off_diag)*sqrt(F); % 1 - 2*die_var
%param_dde.theta_minI(c,p,i) = 0; %-(0.1-2*A); %-2*A_center;
%end

end

% Initialize parameters for image
x0 =x_approx; %x_th.xo; %x_approx;
param_algo.im_approx_thresh = 0.01;


for i =1:3
    xq{i} = x_approx{i};
    x_th.true_best{i} = xq{i};
     x0{i}(abs(x0{i})< param_algo.im_approx_thresh*max(abs(x0{i}(:)))) = 0 ; % 0.01 % 0.05
       mask_app{i} = (abs(x0{i}) > 0);

    param_im.min_x{i} = zeros(size(x0{i}));
    sol_mask{i} = param_im.min_x{i};
    sol_mask{i} = x0{i}(mask_app{i});
%      param_im.min_x{i}(abs(x0{i})>0) = -1e-1*abs(x0{i}(abs(x0{i})>0));
     param_im.max_x{i} = 5e-1*abs(x0{i}(abs(x0{i})>0));

%  if i ~= 1
%       param_im.max_x{i} = zeros(size(x0{i}));
%      param_im.min_x{i}(:) = 5e-2*min(x0{i}(:));
%      param_im.max_x{i}(:) = 5e-2*max(x0{i}(:));
%       param_im.max_x{i}(mask_app{i}) = 4e-2*abs(x0{i}(mask_app{i}));
% % 
%  end

    param_im.min_x{i}(mask_app{i}) = -5e-1*abs(x0{i}(mask_app{i}));
    epsilon{i} = zeros(size(x0{i}));
    x_th.x_approx{i} = x0{i};
    x_th.eps_true{i} = x_th.true_best{i}-x_th.x_approx{i};
     epsilon{i} =  zeros(size(x_th.eps_true{i}));
     x_th.xo{i} = x0{i};
end

% DDE initialisation: /!\ beware amplitude
vec_z = zeros(S2,na, P, 1);
% vec_z(c,:,:,:) = A;
p = floor(P/2) + 1;
U1 = ((randn(S2, na, P, 4) + 1i*randn(S2, na, P,4))/P)*sqrt(F);

U1(:,:,:,1) = A_center*U1(:,:,:,1);
U1(:,:,:,4) = A_center*U1(:,:,:,4);
U1(:,:,:,2) = vec_z.*U1(:,:,:,2);
U1(:,:,:,3) = vec_z.*U1(:,:,:,3);
% % 
% U1(c,:,:,1) = 1;
% U1(c,:,:,4) = 1;
% U1(c,:,:,2) = 0.01;
% U1(c,:,:,3) = 0.01;

U1(:,:,p,:) = 0;
U1(c,:,p,1) = sqrt(F);
U1(c,:,p,4) = sqrt(F);
U1(c,:,p,2) = 0; %sqrt(F)*off_diag;
U1(c,:,p,3) = 0; %sqrt(F)*off_diag;
%U1(:,:,:,4) = U1(:,:,:,1);


%  U1r = max(min(real(U1), reshape(param_dde.theta_maxR, [S2, 1, P])), reshape(param_dde.theta_minR, [S2, 1, P]));
% U1i = max(min(imag(U1), reshape(param_dde.theta_maxI, [S2, 1, P])), reshape(param_dde.theta_minI, [S2, 1, P]));
% U1 = U1r + 1i*U1i;

for i = 1:4
    U1r(:,:,:,i) = max(min(real(U1(:,:,:,i)), reshape(param_dde.theta_maxR(:,:,i), [S2,1,P,1])), reshape(param_dde.theta_minR(:,:,i), [S2,1,P,1]));
    U1i(:,:,:,i) = max(min(imag(U1(:,:,:,i)), reshape(param_dde.theta_maxI(:,:,i), [S2, 1,P,1])), reshape(param_dde.theta_minI(:,:,i), [S2, 1, P,1]));
    U1(:,:,:,i) = U1r(:,:,:,i) + 1i*U1i(:,:,:,i);
end

% U1 = U_th;
U2 = U1;

%%
% U1 = U_th;
% Gt = []; % unsued parameter (left for legacy purpose, to be cleansed later...)
% scale_t = [];
% 
% [~,D1] = computeD_reg(U1, T, Gt, scale_t, T);
% D2 = D1;
% 
%  parfor t = 1:T
%         %     G{t} = createGnufft_T2(D1(:, :, t), D2(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t}); % incorporate misisng antenna positions (evolves w.r.t. t)
%         [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v(:,t), u(:,t)], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
%     end
%     % G0 = cell2mat(G);
%     G0_{1} = cell2mat(G1_');
%     G0_{2} = cell2mat(G2_');
%     G0_{3} = cell2mat(G3_');
%     G0_{4} = cell2mat(G4_');
%     
%     clear G1_
%     clear G2_
%     clear G3_
%     clear G4_
%      
%     
%     
%     for i = 1:4
%         B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
%         Bt_tmp{i} = @(x) (so_fft2_adj((G0_{i})'*x, N, K, sp_scale));
%     end
   
% B_tmp = 1;
% Bt_tmp = 1;


%%

%% Golden search
gs = 0;
% B_tmp = 1;
% Bt_tmp = 1;
if gs
% DEFINE ALL YOUR PARAMETERS FOR YOUR ALGORITHM HERE
% param = ....
% x0 =  ...



% DEFINE YOUR ALGORITHM AS A FUNCTION
% par is the parameter you want to optimize
% algoGS =@(par) - your_algo_GS(x0,param, par);

algoGS = @(par) - joint_imaging_dde_blckCoordinate6_4_gs(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U_th, ...
    U_th, sp_scale, V, param_im, param_dde, param_algo, W, name2, mu,param_pd,L,Lt, B_tmp,Bt_tmp, par, U_th);





% -------------------------------------------------
% golden_search will minimize the first value sent by algoGS 
% -------------------------------------------------
% - if you want to maximize the SNR : 
%   algoGS =@(par) - your_algo_GS(x0,param, par);
%   the first variable sent by your_algo_GS must be the SNR
%   you have to put 'minus' since you want to MAXIMIZE the SNR
%   and golden_search will minimize the value sent by algoGS
% -------------------------------------------------
% - if you want to minimize the MSE : 
%   algoGS =@(par) your_algo_GS(x0,param, par);
%   the first variable sent by your_algo_GS must be the MSE
% -------------------------------------------------

nu_inf = 1e-1 ; % Lowerbound for GS
nu_sup = 1e+3 ; % Upper bound for GS
[nu_opt, save_nu, save_SNR] = golden_search(algoGS, nu_inf, nu_sup) ;
nu = nu_opt ; % your optimized parameter

figure
subplot 211
plot(save_nu, 'o')
xlabel('iterations of Golden search algo')
ylabel('value of the parameter')
subplot 212
plot(save_SNR, 'x')
xlabel('iterations of Golden search algo')
ylabel('value of the corresponding SNR')
end
%%
% DDE estimation (set parameters & estimate DDE + image)

% %diary(sprintf('test_S_%d_method_%d_A_%d_A_center_%d',S,param_pd.dual_fb,A,A_center))

B_tmp = 1;
Bt_tmp = 1;

if strcmp(im_choice,'hydra')
   if param_pd.dual_fb
    param_l1.eta_o = [70, 20, 20];
else
param_l1.eta_o =  [70, 20, 20];
   end
else
if param_pd.dual_fb
    param_l1.eta_o = [500, 400, 400]; %[950, 800, 500]; %[500,400,400];
else
param_l1.eta_o =  [30,20,20]; %[500, 400, 400]; 
end
end


eta__ = param_l1.eta_o;

% eta__ =  [500,400,400]; %eta_o; %[30,20,20]; %[800,400,400]; %[500,400,600]; %eta_o; %[110,100,100]; %,[1120,65,187];
param_algo.eta = eta__;
[U1, U2, epsilon, objective, time_dde, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate6_4(y, Y, x0, ...
    epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, ...
    U2, sp_scale, V, param_im, param_dde, param_algo, W, name2, mu,param_pd,L,Lt, B_tmp,Bt_tmp,eta__,U_th);

%%
for i =1:3
    x{i} = x0{i} + epsilon{i};
 x_reg{i} = x{i}/max(x{i}(:));
 scale_x(i) = sum(x_reg{i}(:).*x_th.x{i}(:))/sum(x_reg{i}(:).^2);
 snr_reg(i) = SNR(scale_x(i)*x_reg{i}, x_th.x{i});
end
% disp(['SNR(x_reg) = ', num2str(snr_reg)]);
% 
% % figure; plot(objective(1:find(objective, 1, 'last')));
 mkdir('results')
 save(['results/joint_imaging_dde_reg_', name2,'.mat'], 'U1', 'U2', 'epsilon', ...
      'x0', 'objective', 'param_dde', 'time_dde', 'error_dirty', 'snr_dirty', 'snr_x', ...
      'snr_reg', 'x_reg', 'snr_approx_die', 'x_approx_die','param_algo','param_pd','L');

  
  
 %% To run primal-dual algorithm for imaging with polarization constraint
 
 
%  %% Load files
% A_center = 5e-2; %5e-2;
% A = 5e-3;
% test_number = A_center;
% off_diag = 0.01;
% rng(2);
% %  util_create_pool_bis(4, []);
% 
% S = 3 %3 %5;
% S_true = S;
% 
% % Data loading
% % name = [cov_type, '_', im_choice, '_', test_number,'_',num2str(S),'_ data_type',num2str(data_type)] %,'DDEs_noise'];
% name = [cov_type, '_', im_choice, '_', num2str(test_number),'_S=',num2str(S),'_ data_type',num2str(data_type)]; %,'_off_diag',num2str(off_diag)];
% 
% name2 = [name, test_aux];
% load(['synth_data_reg_eusipco_', name, '.mat']) % for temporal regularization
% % load(['synth_data_reg_cell_',name2,'.mat']) % without temporal regularization
% 
% load(['joint_imaging_dde_reg_', name2,'.mat']) 
%   load(['imaging_norm_die_', name2, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die')

 %%
 [~,D1] = computeD_reg(U1, T, [], [], T); % estimated DDEs
 [~,D2] = computeD_reg(U2, T, [], [], T); % estimated DDEs
 
 parfor t = 1:T
     %     G{t} = createGnufft_T2(D1(:, :, t), D2(:, :, t), [v(:, t), u(:, t)], K, S, J, W{t}); % incorporate misisng antenna positions (evolves w.r.t. t)
     [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v(:,t), u(:,t)], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
 end
 % G0 = cell2mat(G);
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
 if strcmp(im_choice,'hydra')
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
        %     eta{i} = param_l1.eta_o(i)*weights{i}*1.9/Ax;
    end
end

[x_final, snr] = update_stokes(y, Jeps, tol_x, nIter, A_est, At_est,L, Lt, M, param_l1.eta_o, param_l1,param_pd,x_th, param_dual, Ax, x, eta_);



 
save(['results/joint_imaging_dde_reg_', name2,'_post_processing.mat'], 'x_final','eta_o');


snr_x1 = SNR(x_final{1}, x_th.x{1});
snr_x2 = SNR(x_final{2}, x_th.x{2});
snr_x3 = SNR(x_final{3}, x_th.x{3});

fprintf('%2s\t%2s\t%2s  \n', 'Final SNR_I', 'SNR_Q', 'SNR_U')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e\t%.2e\t%.2e \n',snr_x1, snr_x2, snr_x3);
fprintf('================================================================================\n');
end

