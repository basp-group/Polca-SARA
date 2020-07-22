%% Script to define the necessary parameters for imaging

%%
% -------------------------------------------------------------------------
% Parameters for primal-dual to compute the proximity operator
% -------------------------------------------------------------------------

param_pd.rel_obj = 1e-4; % stopping criterion
param_pd.max_iter = 200; % max number of iterations
param_pd.nnls_init=0;
[param_pd.Ny, param_pd.Nx] = size(x_th.x{1});
param_pd.method = 3;
param_pd.nu1 = 1; % bound on the norm of the operator Psi
% param.nu3 = Ax % bound on the norm of the measurement operator 
param_pd.nu2 = 1; %for ephigraphical projection
param_pd.pol_thresh = 1e-2*param_pd.Ny*param_pd.Nx; % Threshold value for num of pixels not satisfying pol. const.
param_pd.pol_tol = 0; % Tolerance on the polarization constraint: \sqrt(Q^2+U^2+V^2) - I <= tol
param_pd.verbose = 4; % print log or not
param_pd.dual_fb = m_ch; 
param_pd.real = 1;

n_tmp = ['_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp];

%%
% -------------------------------------------------------------------------
% Parameters for regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;

if param_pd.dual_fb
    param_algo.opt_dict = 'SARA'; % 'DB' % 'SARA' %'Dirac' %'DB8'
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
%param_l1.eta_o = [1e-4, 1e-3, 1e-3];
param_l1.approx = 0; %1: If box bounds on the image to be considered

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

 util_create_pool(12);
 
 
Jeps = 800;
tol_x = 1e-4;
nIter = 1;

if (strcmp(im_choice,'hydra') || strcmp(im_choice,'hydra5'))
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

eta_o = param_l1.eta_o;

 
 %%
 % Parameters related to measurements
 
yy = y(:);
s_col = size(yy,1)/4;
y = reshape(yy,[s_col,4]);

param_pd.M = size(y,1); % Number of measurements for each correlator
M = param_pd.M;


if L_ch == 1
    % For linear feeds
    L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
else
    % For circular feeds
    L = [1,0,0,1;0,1,1,0;0,1i,-1i,0];
end
Lt = 0.5*L'; %Adjoint conversion matrix

 