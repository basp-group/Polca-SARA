%% Script to define the necessary parameters for joint calibration and imaging

% -------------------------------------------------------------------------
% regularisation parameters
% -------------------------------------------------------------------------
param_algo.eta = eta_o; % Regularization parameters for the images
param_algo.Jeps = 100; % Number of FB iterations for imaging
param_algo.tol_x = 1e-5; % Stopping criterion for imaging
param_algo.tol_crit = 2e-2; % Global stopping criterion
% parameters for DDE estimation
param_dde.max_it = 20; % Global number of iterations
param_dde.JUtot = 2; % Number of iterations to complete a cycle
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
% param.nu3 = Ax % bound on the norm of the measurement operator
param_pd.nu2 = 1; %for ephigraphical projection
param_pd.nu3 = 1;
param_pd.pol_thresh = 1e-2*param_pd.Ny*param_pd.Nx; % Threshold value for num of pixels not satisfying pol. const.
param_pd.verbose = 4; % print log or not
param_pd.real = 1;


% -------------------------------------------------------------------------
% regularisation (wavelet dictionary)
% -------------------------------------------------------------------------
param_algo.nlevel = 4;

if param_pd.dual_fb
    param_algo.opt_dict ='SARA'; % 'DB' % 'SARA' %'Dirac' %'DB8'
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


eta__ = param_l1.eta_o;

param_algo.eta = eta__;

