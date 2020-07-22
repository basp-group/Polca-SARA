% Code to analyze the results obtained by joint calibration and imaging
% algorithm
function[snr_die2] = fun_results_analyze_die_cal(im_ch, m_ch, L_ch)

addpath data
addpath utils
addpath utils/nufft
addpath utils/alg
addpath lib
addpath results
addpath results_die
addpath ../../Tools
addpath dde_tmp_reg

T = 200;
if im_ch == 1
im_choice = 'cyg_a';
elseif im_ch == 2
im_choice = 'hydra';
else
im_choice = 'hydra5';
end

cov_type = 'vlaa';
na = 27;
S = 5;
S_true = S;
data_type = 3; %1:cal_transfer; 2: cal_transfer+DDEs for off diagonal term; 3: no cal_transfer
A_center = 5e-2; 
A = 5e-2;
A_off = 5e-4;
test_number = 011;
if strcmp(im_choice,'hydra5')
test_number = 511;
end

off_diag = 0;
param_pd.dual_fb = m_ch; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox
S_die = 1;

if L_ch == 1
% For linear feeds
L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
%
else
% For circular feeds
L = [1,0,0,1;0,1,1,0;0,1i,-1i,0];
end

Lt = 0.5*L'; %0.5*(conj(L))'; %Adjoint conversion matrix


snr_die2 = zeros(5,3);

for seed = 1:5
rng(seed);

% Loading data and results
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


% load(['/data/jbirdi/Pol_cal/final_cirrus/data/synth_data_reg_eusipco_', name,'.mat'])
load(['data/synth_data_reg_cell_', name, '.mat'])


yy = y(:);
s_col = size(yy,1)/4; 
y = reshape(yy,[s_col,4]);
M = 70200; %param_pd.M;
N = size(x_th.x{1}); %size(x_th.x{1},1) * size(x_th.x{1},2);

load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die')

n_tmp_ = ['_new_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];
load(['results/joint_imaging_die_reg_', name2,'.mat'])

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

for i =1:3
    x{i} = x0{i} + epsilon{i};
    x_reg{i} = x{i}/max(x{i}(:));
    scale_x(i) = sum(x_reg{i}(:).*x_th.x{i}(:))/sum(x_reg{i}(:).^2);
    snr_die_cal(i) = SNR(scale_x(i)*x_reg{i}, x_th.x{i});
    x_die_cal{seed,i} = x{i};
end

% IMAGES
% True images: x_th.x
% Images obtained from calibrated DIEs: x
% Scaled images obtained from calibrated DIEs: x_reg


% CALIBRATION 
% True DDEs: U_th
% Calibrated DIEs: D1, D2

snr_die2(seed,:) = snr_die_cal;


% Weighted l2 norm of residual images
% With calibrated DIEs

%Create convolution matrix G

[~,D_true] = computeD_reg(U_th, T, [], [], T);

    for t = 1:T
        [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v{t}, u{t}], K, S_die, J, W{t});
        [G1{t},G2{t},G3{t},G4{t}] = createGnufft_T2_parallel_T(D_true{t},D_true{t},[v{t}, u{t}], K, S, J, W{t});
    end

    % G0 = cell2mat(G);
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
        B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
        Bt_tmp{i} = @(x) (so_fft2_adj((G0_{i})'*x, N, K, sp_scale));

        A_true{i} = @(x) G0{i}*so_fft2(x, K, sp_scale);
        At_true{i} = @(x) (so_fft2_adj((G0{i})'*x, N, K, sp_scale));
    end

    lambda_scale = 1;
    
    Ax(seed) = 2*lambda_scale*op_norm_stokes_RIME(B_tmp,Bt_tmp, size(x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
    Ax_true(seed) = 2*lambda_scale*op_norm_stokes_RIME(A_true,At_true, size(x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case

    % For images obtained from calibration algo: x
    b_x = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
    
    [~,y_x] = Phi(B_tmp, b_x, M);
    [~,res_x1] = Phit(Bt_tmp,y_x-y,1,L,Lt);
    
    [~,y_xt] = Phi(A_true, b_x, M);
    [~,res_x1t] = Phit(At_true,y_xt-y,1,L,Lt);
    
   
    
    for i = 1:3
        res_x{seed,i} = res_x1{i};
        res_xtrue{seed,i} = res_x1t{i};
        
        % Weighted l2 norm of residual images
        w_resx(seed,i) = norm(res_x{seed,i},2)/sqrt(prod(N)); % with estimated DDEs
        w_resxt(seed,i) = norm(res_xtrue{seed,i},2)/sqrt(prod(N)); % with true DDEs
    end
    
    % Dynamic range calculation
    for i =1:3
        dr(seed,i) = (Ax(seed)*max(x{i}(:)))/w_resx(seed,i);
        
        dr_t(seed,i) = (Ax_true(seed)*max(x{i}(:)))/w_resxt(seed,i);
        
    end

end


% Compute the mean and the best values of the metrics over 5 simulations

mean_snr_die_cal = mean(snr_die2);


mean_dr = mean(dr);
mean_dr_t = mean(dr_t);

mean_w_resx = mean(w_resx);
mean_w_resxt = mean(w_resxt);

   
[max_snr_die2, max_snr_die2_seed] = max(snr_die2);


[max_dr, max_dr_seed] = max(dr);
[max_dr_t, max_dr_t_seed] = max(dr_t);


[min_w_resx, min_w_resx_seed] = min(w_resx);
[min_w_resxt, min_w_resxt_seed] = min(w_resxt);
                                       
                                       %%
                                       % Save images
if im_ch == 1
seed = 2;
else
seed = 1;
end
                                       
                                       if param_pd.dual_fb == 0
                                       name_save = sprintf('reg_die_%s_%s.mat',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       save(['results_images/die_cal/saved_files/',name_save],'L','snr_die2','x_die_cal','res_x','res_xtrue',...
                                            'w_resx','w_resxt','dr','dr_t','mean_snr_die_cal','mean_dr','mean_dr_t',...
                                            'mean_w_resx','mean_w_resxt','Ax','Ax_true','seed');
                                       
                                       
                                       name_fig1 = sprintf('rec_I_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps1 = sprintf('rec_I_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       name_fig2 = sprintf('rec_Q_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps2 = sprintf('rec_Q_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       name_fig3 = sprintf('rec_U_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps3 = sprintf('rec_U_reg_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       else
                                       name_save = sprintf('const_die_%s_%s.mat',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       save(['results_images/die_cal/saved_files/',name_save],'L','snr_die2','x_die_cal','res_x','res_xtrue',...
                                            'w_resx','w_resxt','dr','dr_t','mean_snr_die_cal','mean_dr','mean_dr_t',...
                                            'mean_w_resx','mean_w_resxt','Ax','Ax_true','seed');
                                       
                                       
                                       name_fig1 = sprintf('rec_I_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps1 = sprintf('rec_I_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       name_fig2 = sprintf('rec_Q_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps2 = sprintf('rec_Q_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       
                                       name_fig3 = sprintf('rec_U_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       name_eps3 = sprintf('rec_U_const_die_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
                                       end
                                       
                                       i = 1;
                                       
                                       
                                       if strcmp(im_choice,'hydra')
                                       cmin = -3.4;
                                       else
                                       cmin = -3.2;
                                       end
                                       %     cmin = -3.2; %-3.4% min(log10(abs(x_th.x{i}(:))));
                                       cmax = max(log10(abs(x_th.x{i}(:))));
                                       
                                       
                                       
                                       figure, imagesc(log10(abs(x_die_cal{seed,i}))),colormap('jet'), caxis([cmin cmax]), axis image, axis off; set(gca, 'Position',[0 0 1 0.99],'visible','off')
                                       savefig(sprintf('results_images/die_cal/save_fig/%s.fig',name_fig1));
                                       saveas(gcf,sprintf('results_images/die_cal/save_eps/%s.eps',name_eps1),'epsc');
                                       
                                       close all
                                       
                                       i = 2;
                                       
                                       if strcmp(im_choice,'hydra')
                                       cmin = -4.4;
                                       else
                                       cmin = -3.5;
                                       end
                                       
                                       
                                       figure, imagesc(log10(abs(x_die_cal{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
                                       savefig(sprintf('results_images/die_cal/save_fig/%s.fig',name_fig2));
                                       saveas(gcf,sprintf('results_images/die_cal/save_eps/%s.eps',name_eps2),'epsc');
                                       
                                       close all;
                                       
                                       i = 3;
                                       %cmin = -4.4; %-3.5; %min(log10(abs(x_th.x{i}(:))));
                                       cmax = max(log10(abs(x_th.x{i}(:))));
                                       
                                       figure, imagesc(log10(abs(x_die_cal{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
                                       savefig(sprintf('results_images/die_cal/save_fig/%s.fig',name_fig3));
                                       saveas(gcf,sprintf('results_images/die_cal/save_eps/%s.eps',name_eps3),'epsc');
                                       
                                       close all;
                                       
                                       fprintf('%2s  \n', 'Mean SNR_die_cal')
                                       fprintf('--------------------------------------------------------------------------------\n');
                                       fprintf('%.2e \n',mean_snr_die_cal);
                                       fprintf('================================================================================\n');
                                       

                                       fprintf('%2s  \n', 'Mean DR_die_cal')
                                       fprintf('--------------------------------------------------------------------------------\n');
                                       fprintf('%.2e \n',mean_dr);
                                       fprintf('================================================================================\n');

                                       

                                       
end

