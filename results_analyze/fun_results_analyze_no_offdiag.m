% Code to analyze the results obtained by joint calibration and imaging
% algorithm
function[snr_die,snr_final] = fun_results_analyze_no_offdiag(im_ch,m_ch,L_ch)
%close all; clear all;

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

na = 27;
cov_type = 'vlaa';
S = 5;
S_true = S;
data_type = 3; %1:cal_transfer; 2: cal_transfer+DDEs for off diagonal term; 3: no cal_transfer
A_center = 5e-2; 
A = 5e-2;
A_off = 5e-4;
test_number = 411; %011; 
off_diag = 0;
param_pd.dual_fb = m_ch; %0: run primal dual to compute prox; 1: run dual forward backward to compute prox


if strcmp(im_choice,'hydra5')
test_number = 511;
end


if L_ch == 1
% For linear feeds
  L = [1,0,0,1;1,0,0,-1;0,1,1,0]; % Conversion matrix
%
else
% For circular feeds
 L = [1,0,0,1;0,1,1,0;0,1i,-1i,0];
end
Lt = 0.5*L'; %0.5*(conj(L))'; %Adjoint conversion matrix


snr_die = zeros(5,3);
snr_cal = zeros(5,3);
snr_final = zeros(5,3);

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


 load(['synth_data_reg_eusipco_', name,'.mat'],'D_th','F','J','K','N','Omega','P','U_th','V','W','Y','sp_scale','u','v','x_th','y');
%load(['data/synth_data_reg_eusipco_', name,'.mat'])

yy = y(:);
s_col = size(yy,1)/4; 
y = reshape(yy,[s_col,4]);
N = size(x_th.x{1});

load(['imaging_norm_die_', name_, '.mat'], 'eta_o', 'x_approx','snr_approx_die_unscaled','xq','x_approx_die','snr_approx_die')

if param_pd.dual_fb == 0
n_tmp_ = ['_no_off_diag_new_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];
load(['results/joint_imaging_dde_reg_', name2,'.mat'])
else
 n_tmp_ = ['_no_offdiag_new_with_const_method_',num2str(param_pd.dual_fb)];
name2 = [name, n_tmp_];
%name2 = [name2_, n_tmp];
load(['results/joint_imaging_dde_reg_mod_pdfb_', name2,'.mat'])
end

load(['results/joint_imaging_dde_reg_', name2,'_post_processing.mat'])

M = param_pd.M;

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

for i =1:3
    x{i} = x0{i} + epsilon{i};
    x_ff{i} = x_final{i}/max(x_final{i}(:));
    scale_ff(i) = sum(x_ff{i}(:).*x_th.x{i}(:))/sum(x_ff{i}(:).^2);
    snr_scale_final(i) = SNR(scale_ff(i)*x_ff{i}, x_th.x{i});
end

% IMAGES
% True images: x_th.x
% Images obtained from normalized DIEs: x_approx
% Scaled images obtained from normalized DIEs: x_approx_die
% Images obtained from calibration algo: x
% Scaled images obtained from calibration algo: x_reg
% Images obtained from final imaging step: x_final
% Scaled images obtained from final imaging step: x_ff

% CALIBRATION 
% True DDEs: U_th
% Calibrated DDEs: U1, U2

snr_die(seed,:) = snr_approx_die;
snr_cal(seed,:)= snr_reg;
snr_final(seed,:) = snr_scale_final;

for i =1:3
    x_die{seed,i} = x_approx{i};
    x_cal{seed,i} = x{i};
    x_f{seed,i} = x_final{i};
end



% Weighted l2 norm of residual images
% With calibrated DDEs

%Create convolution matrix G
[~,D1] = computeD_reg(U1, T, [], [], T);
[~,D2] = computeD_reg(U2, T, [], [], T);

[~,D_true] = computeD_reg(U_th, T, [], [], T);

for t = 1:T
    D_approx{t} = ones(4,1,na);
    D_approx{t}(2,1,:) = 0; %(0.01 + 2*A); %0;
    D_approx{t}(3,1,:) = 0; %(0.01 + 2*A); % 0;
end


parfor t = 1:T
    [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(D1{t}, D2{t}, [v(:,t), u(:,t)], K, S, J, W{t});
    [G1{t},G2{t},G3{t},G4{t}] = createGnufft_T2_parallel_T(D_true{t},D_true{t},[v(:,t), u(:,t)], K, S, J, W{t});
    [G1__{t},G2__{t},G3__{t},G4__{t}] = createGnufft_T2_parallel_T(D_approx{t}, D_approx{t}, [v(:,t), u(:,t)], K, 1, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format

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
    
    G0__{1} = cell2mat(G1__');
    G0__{2} = cell2mat(G2__');
    G0__{3} = cell2mat(G3__');
    G0__{4} = cell2mat(G4__');

    clear G1_; clear G1; clear G1__
    clear G2_; clear G2; clear G2__
    clear G3_; clear G3; clear G3__
    clear G4_; clear G4; clear G4__



    for i = 1:4
        B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, sp_scale);
        Bt_tmp{i} = @(x) (so_fft2_adj((G0_{i})'*x, N, K, sp_scale));

        A_true{i} = @(x) G0{i}*so_fft2(x, K, sp_scale);
        At_true{i} = @(x) (so_fft2_adj((G0{i})'*x, N, K, sp_scale));
        
        A_die{i} = @(x) G0__{i}*so_fft2(x, K, sp_scale);
        At_die{i} = @(x) (so_fft2_adj((G0__{i})'*x, N, K, sp_scale)); % for circular feed - no real

    end

    lambda_scale = 1;
    
    Ax(seed) = 2*lambda_scale*op_norm_stokes_RIME(B_tmp,Bt_tmp, size(x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
    Ax_true(seed) = 2*lambda_scale*op_norm_stokes_RIME(A_true,At_true, size(x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case
    Ax_die(seed) = 2*lambda_scale*op_norm_stokes_RIME(A_die,At_die, size(x{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case

    % For images obtained from normalized DIEs: x_approx
    b_xa = conv_stokes_to_bright([x_approx{1}(:),x_approx{2}(:),x_approx{3}(:)],L,size(x_approx{1}));
    
    [~,y_xa] = Phi(A_die, b_xa, M);
    [~,res_x1a] = Phit(At_die,y_xa-y,1,L,Lt);
    
    [~,y_xat] = Phi(A_true, b_xa, M);
    [~,res_xat] = Phit(At_true,y_xat-y,1,L,Lt);
    
    % For images obtained from calibration algo: x
    b_x = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
    
    [~,y_x] = Phi(B_tmp, b_x, M);
    [~,res_x1] = Phit(Bt_tmp,y_x-y,1,L,Lt);
    
    [~,y_xt] = Phi(A_true, b_x, M);
    [~,res_x1t] = Phit(At_true,y_xt-y,1,L,Lt);
    
    % For images obtained from final imaging step: x_final
    b_xfinal = conv_stokes_to_bright([x_final{1}(:),x_final{2}(:),x_final{3}(:)],L,size(x_final{1}));
    
    [~,y_xfinal] = Phi(B_tmp, b_xfinal, M);
    [~,res_x2] = Phit(Bt_tmp,y_xfinal-y,1,L,Lt);
    
    [~,y_xfinalt] = Phi(A_true, b_xfinal, M);
    [~,res_x2t] = Phit(At_true,y_xfinalt-y,1,L,Lt);
    
    for i = 1:3
        res_xa{seed,i} = res_x1a{i};
        res_xatrue{seed,i} = res_xat{i};
        
        res_x{seed,i} = res_x1{i};
        res_xf{seed,i} = res_x2{i};
        
        res_xtrue{seed,i} = res_x1t{i};
        res_xftrue{seed,i} = res_x2t{i};
        
        % Weighted l2 norm of residual images
        w_resx(seed,i) = norm(res_x{seed,i},2)/sqrt(prod(N)); % with estimated DDEs
        w_resxf(seed,i) = norm(res_xf{seed,i},2)/sqrt(prod(N));
        w_resxa(seed,i) = norm(res_xa{seed,i},2)/sqrt(prod(N));

        w_resxt(seed,i) = norm(res_xtrue{seed,i},2)/sqrt(prod(N)); % with true DDEs
        w_resxft(seed,i) = norm(res_xftrue{seed,i},2)/sqrt(prod(N));
        w_resxat(seed,i) = norm(res_xatrue{seed,i},2)/sqrt(prod(N));
        
        
    end
    
    % Dynamic range calculation
    for i =1:3
        dr(seed,i) = (Ax(seed)*max(x{i}(:)))/w_resx(seed,i);
        dr_f(seed,i) = (Ax(seed)*max(x_final{i}(:)))/w_resxf(seed,i);
        dr_a(seed,i) = (Ax_die(seed)*max(x_approx{i}(:)))/w_resxa(seed,i);

        dr_t(seed,i) = (Ax_true(seed)*max(x{i}(:)))/w_resxt(seed,i);
        dr_ft(seed,i) = (Ax_true(seed)*max(x_final{i}(:)))/w_resxft(seed,i);
        dr_at(seed,i) = (Ax_true(seed)*max(x_approx{i}(:)))/w_resxat(seed,i);

    end

end


% Compute the mean and the best values of the metrics over 5 simulations


mean_snr_die = mean(snr_die);
mean_snr_cal = mean(snr_cal);
mean_snr_final = mean(snr_final);

mean_dr = mean(dr);
mean_dr_f = mean(dr_f);
mean_dr_a = mean(dr_a);
mean_dr_t = mean(dr_t);
mean_dr_ft = mean(dr_ft);
mean_dr_at = mean(dr_at);

mean_w_resx = mean(w_resx);
mean_w_resxf = mean(w_resxf);
mean_w_resxa = mean(w_resxa);
mean_w_resxt = mean(w_resxt);
mean_w_resxft = mean(w_resxft);
mean_w_resxat = mean(w_resxat);

   
[max_snr_die, max_snr_die_seed] = max(snr_die);
[max_snr_cal, max_snr_cal_seed] = max(snr_cal);
[max_snr_final, max_snr_final_seed] = max(snr_final);


[max_dr, max_dr_seed] = max(dr);
[max_dr_f, max_dr_f_seed] = max(dr_f);
[max_dr_a, max_dr_a_seed] = max(dr_a);

[max_dr_t, max_dr_t_seed] = max(dr_t);
[max_dr_ft, max_dr_ft_seed] = max(dr_ft);
[max_dr_at, max_dr_at_seed] = max(dr_at);


[min_w_resx, min_w_resx_seed] = min(w_resx);
[min_w_resxf, min_w_resxf_seed] = min(w_resxf);
[min_w_resxa, min_w_resxa_seed] = min(w_resxa);

[min_w_resxt, min_w_resxt_seed] = min(w_resxt);
[min_w_resxft, min_w_resxft_seed] = min(w_resxft);
[min_w_resxat, min_w_resxat_seed] = min(w_resxat);


%%
% Save images
if im_ch == 1
seed = 2;
else
seed = 1;
end
                                      
if param_pd.dual_fb == 0
    name_save = sprintf('reg_no_offdiag_%s_%s.mat',num2str(im_choice),num2str(A_center),num2str(test_number));

save(['results_images/saved_files/',name_save],'L','snr_die','snr_cal','snr_final','x_die','x_cal','x_f','res_x','res_xf','res_xa','res_xtrue','res_xftrue','res_xatrue',...
    'w_resx','w_resxf','w_resxa','w_resxt','w_resxft','w_resxat','dr','dr_f','dr_a','dr_t',...
'dr_ft','dr_at','mean_snr_die','mean_snr_cal','mean_snr_final','mean_dr','mean_dr_f','mean_dr_t','mean_dr_a','mean_dr_ft','mean_dr_at',...
'mean_w_resx','mean_w_resxf','mean_w_resxa','mean_w_resxt','mean_w_resxft','mean_w_resxat','Ax','Ax_true','Ax_die','seed');

    
    name_fig1 = sprintf('rec_I_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps1 = sprintf('rec_I_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    
    name_fig2 = sprintf('rec_Q_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps2 = sprintf('rec_Q_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    
    name_fig3 = sprintf('rec_U_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps3 = sprintf('rec_U_reg_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    
else
        name_save = sprintf('const_no_offdiag_%s_%s.mat',num2str(im_choice),num2str(A_center),num2str(test_number));

 save(['results_images/saved_files/',name_save],'L','snr_die','snr_cal','snr_final','x_die','x_cal','x_f','res_x','res_xf','res_xa','res_xtrue','res_xftrue','res_xatrue',...
    'w_resx','w_resxf','w_resxa','w_resxt','w_resxft','w_resxat','dr','dr_f','dr_a','dr_t',...
'dr_ft','dr_at','mean_snr_die','mean_snr_cal','mean_snr_final','mean_dr','mean_dr_f','mean_dr_t','mean_dr_a','mean_dr_ft','mean_dr_at',...
'mean_w_resx','mean_w_resxf','mean_w_resxa','mean_w_resxt','mean_w_resxft','mean_w_resxat','Ax','Ax_true','Ax_die','seed');

 
    name_fig1 = sprintf('rec_I_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps1 = sprintf('rec_I_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    
    name_fig2 = sprintf('rec_Q_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps2 = sprintf('rec_Q_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    
    name_fig3 = sprintf('rec_U_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
    name_eps3 = sprintf('rec_U_const_no_offdiag_%s_%s',num2str(im_choice),num2str(A_center),num2str(test_number));
end

i = 1;


if strcmp(im_choice,'hydra')
    cmin = -3.4;
else
    cmin = -3.2;
end
%     cmin = -3.2; %-3.4% min(log10(abs(x_th.x{i}(:))));
    cmax = max(log10(abs(x_th.x{i}(:))));
    figure, imagesc(log10(abs(x_th.x{i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/true_I_%s.fig',num2str(im_choice)));
    saveas(gcf,sprintf('results_images/save_eps/true_I_%s.eps',num2str(im_choice)),'epsc');
    
    figure, imagesc(log10(abs(x_die{seed,i}))),colormap('jet'), caxis([cmin cmax]), axis image, axis off; set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_die.fig',name_fig1));
    saveas(gcf,sprintf('results_images/save_eps/%s_die.eps',name_eps1),'epsc');
    
    figure, imagesc(log10(abs(x_cal{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_cal.fig',name_fig1));
    saveas(gcf,sprintf('results_images/save_eps/%s_cal.eps',name_eps1),'epsc');
    
    figure, imagesc(log10(abs(x_f{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_final.fig',name_fig1));
    saveas(gcf,sprintf('results_images/save_eps/%s_final.eps',name_eps1),'epsc');
    
    close all
    
    i = 2;
    
if strcmp(im_choice,'hydra')
    cmin = -4.4;
else
    cmin = -3.5;
end
%     cmin = -3.5; %-4.4; %min(log10(abs(x_th.x{i}(:))));
    cmax = max(log10(abs(x_th.x{i}(:))));
    figure, imagesc(log10(abs(x_th.x{i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/true_Q_%s.fig',num2str(im_choice)));
    saveas(gcf,sprintf('results_images/save_eps/true_Q_%s.eps',num2str(im_choice)),'epsc');
    
    figure, imagesc(log10(abs(x_die{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_die.fig',name_fig2));
    saveas(gcf,sprintf('results_images/save_eps/%s_die.eps',name_eps2),'epsc');
    
    figure, imagesc(log10(abs(x_cal{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_cal.fig',name_fig2));
    saveas(gcf,sprintf('results_images/save_eps/%s_cal.eps',name_eps2),'epsc');
    
    figure, imagesc(log10(abs(x_f{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_final.fig',name_fig2));
    saveas(gcf,sprintf('results_images/save_eps/%s_final.eps',name_eps2),'epsc');
    close all;
    
    i = 3;
    %cmin = -4.4; %-3.5; %min(log10(abs(x_th.x{i}(:))));
    cmax = max(log10(abs(x_th.x{i}(:))));
    figure, imagesc(log10(abs(x_th.x{i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/true_U_%s.fig',num2str(im_choice)));
    saveas(gcf,sprintf('results_images/save_eps/true_U_%s.eps',num2str(im_choice)),'epsc');
    
    figure, imagesc(log10(abs(x_die{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_die.fig',name_fig3));
    saveas(gcf,sprintf('results_images/save_eps/%s_die.eps',name_eps3),'epsc');
    
    figure, imagesc(log10(abs(x_cal{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_cal.fig',name_fig3));
    saveas(gcf,sprintf('results_images/save_eps/%s_cal.eps',name_eps3),'epsc');
    
    figure, imagesc(log10(abs(x_f{seed,i}))), axis image, axis off, colormap('jet'), caxis([cmin cmax]); set(gca, 'Position',[0 0 1 0.99],'visible','off')
    savefig(sprintf('results_images/save_fig/%s_final.fig',name_fig3));
    saveas(gcf,sprintf('results_images/save_eps/%s_final.eps',name_eps3),'epsc');
    close all;
    



fprintf('%2s  \n', 'Mean SNR_die')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_snr_die);
fprintf('================================================================================\n');

fprintf('%2s  \n', 'Mean SNR_cal')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_snr_cal);
fprintf('================================================================================\n');

fprintf('%2s  \n', 'Mean SNR_final')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_snr_final);
fprintf('================================================================================\n');

fprintf('%2s  \n', 'Mean DR_die')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_dr_a);
fprintf('================================================================================\n');

fprintf('%2s  \n', 'Mean DR_cal')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_dr);
fprintf('================================================================================\n');

fprintf('%2s  \n', 'Mean DR_final')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%.2e \n',mean_dr_f);
fprintf('================================================================================\n');


    end


