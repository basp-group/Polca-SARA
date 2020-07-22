% function [St_im, L1_v, L1_vp, L2_vv, L2_vp, rel_sol_norm_change, sol_v, snr_v, no_sub_itr_v, dual_var1, dual_var3, sol_best_bound_v, sol_reweight_v] = pdfb_par_rescaled_precond_wave_par_stokes_rw(y, epsilont, epsilonts, A, At, T, pU, W, Psi, Psit, Psiw, Psitw, Gw,stokes_P, param)
  function [St_im, param_dual, l1] = pdfb_stokes_imaging_cal_ddes(Z, stokes_P, param_l1, param,x_th, param_dual, eta_)
% 
% Z = eps_tmp;
% stokes_P = 3;
% param = param_pd;

clear snr
clear rel_sol_norm_change
%  util_create_pool_bis(4, []);

%% Functions used in proximity operators
negt = @(z) max(real(z), 0); % thresholding negative values
soft = @(z, T) sign(z) .* max(abs(z)-T, 0); %soft thresholding operator
pos_soft = @(z, T) max(real(z)-T, 0); % Positive soft thresholding operator
% param.Ps = 1;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));
method = param.method;
R=1;
Nx = param.Nx;
Ny = param.Ny;
N= Nx*Ny;
Ps = param.Ps;
method = param.method;
evl = param.evl;
Psi = param_l1.Psi;
Psit = param_l1.Psit;
P = 3;
dual_var1 = param_dual.dual_var1;
dual_var2 = param_dual.dual_var2;
dual_tilde2 = param_dual.dual_tilde2;
% 
% Psitw = cell(1,stokes_P);
% Psiw = cell(1,stokes_P);
% 
% for i = 1:stokes_P
%     Psitw{i} = Psit;
%     Psiw{i} = Psi;
% end
% Psitw = param_dual.Psitw;
% Psiw = param_dual.Psiw;
% St_im = Z;
St_im = param_dual.St_im;
% St_im{4} = param_dual.St_im{4};
% St_im{5} = param_dual.St_im{5};

% param_l1.eta_o = [1e-4, 1e-5, 1e-5];

% if (param.method == 0 || param.nnls_init == 0) && n_it == 1
% def_var_real_data;
% 
% %init_var_real_data;
% for i = 1:5      
%   St_im{i} = zeros(Ny,Nx);
% %  dual_tilde1{i} = zeros(size(param_sim_data.Psit(St_im{i})));
% for k=1:Ps
%  dual_var1{i,k} = zeros(size(Psitw{k}(St_im{i})));
% end
%  dual_tilde2{i} = zeros(size(St_im{i}));
%  dual_var2{i} = zeros(size(St_im{i}));
%  dual_tilde3{i} = zeros(size(St_im{i}));
%  dual_var3{i} = zeros(param.M,1); 
% end
% end
%


R = 1;
% % Change the cell formatting

% if param.nnls_init
%     for i=1:stokes_P
%         %         St_im{i} = result_st_nnls.sol{i};
%         for k = 1:param.Ps
%             d_val = abs(Psit{k}(reshape(St_im{i},Ny,Nx)));
%             d_val_max = max(d_val(:));
%             if d_val_max == 0
%                 d_val_max =1;
%             end
%             weights{i,k} = param.reweight_alpha(k) ./ (param.reweight_alpha(k) + (d_val./d_val_max));
%             weights{i,k}(d_val > d_val_max * param.reweight_abs_of_max) = 0;
%             eta{i,k} = param.eta_o(i)*weights{i,k};
%             param.reweight_alpha(k) = param.reweight_alpha_ff(k) .* param.reweight_alpha(k);
%         end
%     end
%     
% else
%     for i = 1:stokes_P
%     for k=1:param.Ps
%          weights{i,k} =  ones(size(Psitw{k}((zeros(Ny,Nx))), 1), 1);
%           eta{i,k} = param_l1.eta_o(i)*weights{i,k};
%     end
%     end
% end

eta = eta_;

 
parfor k = 1:Ps
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end
Psiw = [];
Psiwt = [];
% util_create_pool(Ps);



if ~isfield(param, 'global_stop_bound'), param.global_stop_bound = 1; end
if ~isfield(param, 'rel_obj'), param.rel_obj = 1e-4; end
if ~isfield(param, 'max_iter'), param.max_iter = 200; end
if ~isfield(param, 'nu1')
    param.nu1 = ones(Ps, 1);
else
    if numel(param.nu1) == 1
        param.nu1 = ones(Ps, 1) * param.nu1;
    end
end

if ~isfield(param, 'nu2') % For epigraphical projection
param.nu2 = 1;
end


if ~isfield(param, 'sigma1'), param.sigma1 = 1./param.nu1; end
if ~isfield(param, 'sigma2'), param.sigma2 = 1./param.nu2; end
if ~isfield(param, 'sigma3'), param.sigma3 = 1./param.nu3; end

if ~isfield(param, 'gamma'), param.gamma = 1e-3; end
if ~isfield(param, 'tau'), param.tau = 0.49; end

if method == 1
    param.tau = 1/(0.5 + 1);
else
    param.tau = 1/(0.5 +1 +1);
end

 %% set up log variables
L1_v = zeros(param.max_iter, 1);

%% store useful variables
% step size for the dual variables
sigma1 = param.sigma1;
sigma2 = param.sigma2;
sigma3 = param.sigma3;

% step size primal 
tau = param.tau;

 %%
%****************************************************************************************%
%****************************************************************************************%

% Algorithm

g = zeros(size(St_im{1})); %gradient part 

% Initialize variables 
for i=1:stokes_P
St_im_old{i} = St_im{i};
g1{i} = zeros(size(St_im{i}));
g2{i} = zeros(size(St_im{i}));
g3{i} = zeros(size(St_im{i}));
prev_xsol{i} = 2*St_im{i}-St_im_old{i};
end

% St_im = x_th.im_true;

%%

for t = 1:param.max_iter
tm = tic;
rm(t) = 0;

for i=1:stokes_P

%% Primal update
St_im_old{i} = St_im{i};

if method == 1
    St_im{i} = real((St_im{i} - tau*(reshape(g1{i},Ny,Nx) + (reshape((St_im{i}-Z{i}),Ny,Nx)))));
else
    St_im{i} = real((St_im{i} - tau*(reshape(g1{i},Ny,Nx) + reshape(g2{i},Ny,Nx)+ reshape(St_im{i}-Z{i},Ny,Nx))));
end
% if i == 1
%     St_im{i} = negt(St_im{i});
% end


%%%%%% To constrain the values within a box
if param_l1.approx == 1
 sol_mask{i} = St_im{i}(param_l1.mask_app{i}) ;
if i == 1
    St_im{i}(St_im{i}<param_l1.min_x{i}) = param_l1.min_x{i}(St_im{i}<param_l1.min_x{i});
else
    min_x{i} = param_l1.min_x{i}(param_l1.mask_app{i});
    sol_mask{i}(sol_mask{i}< min_x{i}) = min_x{i}(sol_mask{i}<min_x{i});
end
sol_mask{i}(sol_mask{i}>param_l1.max_x{i}) =param_l1.max_x{i}(sol_mask{i}>param_l1.max_x{i}) ;
St_im{i}(param_l1.mask_app{i}) = sol_mask{i} ;
elseif i == 1     
    St_im{i} = negt(St_im{i});
 end 

if param.real 
    St_im{i} = real(St_im{i}); 
end
%%%%%
prev_xsol{i} = 2*St_im{i}-St_im_old{i};


%% Relative change in the estimates

norm_prevsol(i) = norm(St_im_old{i}(:));
if (norm_prevsol(i) == 0)
    rel_sol_norm_change{i}(t) = 1;
     snr{i}(t) = 0;
     rm(t) = 1;
else
    rel_sol_norm_change{i}(t) = norm(St_im{i}(:) - St_im_old{i}(:))/norm(St_im{i}(:));
    rm(t) = max(rel_sol_norm_change{i}(t),rm(t));
    snr{i}(t) = SNR(St_im{i}, x_th.im_true{i});
end

% end

%% Dual variables update

%% L1 prox update
%%%%%#### Parallelize for the sparsity bases
%     dual_tilde1{i} = dual_var1{i} + sigma1*Psitw{i}(2*St_im{i}-St_im_old{i});
%     dual_var1{i} = dual_tilde1{i} - sigma1*soft(dual_tilde1{i}./sigma1, eta{i}/sigma1);
%     g1{i} = Psiw{i}(dual_var1{i});
end
%

for k = 1:Ps

%       [dual_var1, g1, norm_l1] =  run_par_waverec2(dual_var1{1,k},Psit{k}, Psi{k}, prev_xsol{1}, sigma1(k), St_im{1}, eta{1,k},x_th.psit_xo{1,k});
% end

f1(k) = parfeval(@run_par_waverec2,3,dual_var1{1,k},Psit{k}, Psi{k}, prev_xsol{1}, sigma1(k), St_im{1}, eta{1,k},x_th.psit_xo{1,k});
end

if stokes_P>1
i=2;
for k = 1:Ps
    f2(k) = parfeval(@run_par_waverec2,3,dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{2}, eta{2,k}, x_th.psit_xo{i,k});
end

i=3;
for k = 1:Ps
%         [~] =  run_par_waverec2(dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{i}, eta{i,k},x_th.psit_xo{i,k});

    f3(k) = parfeval(@run_par_waverec2,3,dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{3}, eta{3,k}, x_th.psit_xo{i,k});
end
end

%% Epigraphical projection

if method == 3 % Sparsity + polarization constraint
i = 4;

%Primal update
St_im_old{i} = St_im{i};
St_im_old{5} = St_im{5};
g11{i} = [St_im{i}(:),St_im{i+1}(:)];
g22{i} = tau*sigma2*[dual_var2{i}(:), dual_var2{i+1}(:)];
St{i} = Pv_tol(g11{i}-g22{i}, param.pol_tol);
% St{i} = Pv(g11{i}-g22{i});
St_im{i} = reshape(St{i}(:,1),Ny,Nx);
St_im{i+1} = reshape(St{i}(:,2),Ny,Nx);
prev_xsol{i} = 2*St_im{i}-St_im_old{i};
prev_xsol{5} = 2*St_im{5}-St_im_old{5};

%Dual update
for i=1:stokes_P+2
    dual_tilde2{i} = dual_var2{i} + prev_xsol{i};
end

% par
% vec = [dual_tilde2{1}(:),dual_tilde2{4}(:)];
vec = [dual_tilde2{1}(:)+param_l1.xA{1}(:),dual_tilde2{4}(:)];  % Modified to incorporate x_0 (known bright sources image)
epi_h{1} = vec - Proj_epih1(1,vec);
% epi_h{1}(:,1) = epi_h{1}(:,1) - param_l1.xA{1}(:);

vec = [dual_tilde2{2}(:)+param_l1.xA{2}(:),dual_tilde2{3}(:)+param_l1.xA{3}(:),dual_tilde2{5}(:)]; % Modified to incorporate x_0 (known bright sources image)
epi_h{2} = vec - Proj_epih2(1,vec);
% epi_h{2}(:,1) = epi_h{2}(:,1) - param_l1.xA{2}(:);
% epi_h{2}(:,2) = epi_h{2}(:,2) - param_l1.xA{3}(:);


% Update variables
dual_var2{1} = reshape(epi_h{1}(:,1),Ny,Nx);
dual_var2{4} = reshape(epi_h{1}(:,2),Ny,Nx);

dual_var2{2} = reshape(epi_h{2}(:,1),Ny,Nx);
dual_var2{3} = reshape(epi_h{2}(:,2),Ny,Nx);
dual_var2{5} = reshape(epi_h{2}(:,3),Ny,Nx);

for i=1:stokes_P+2
    g2{i} = sigma2*dual_var2{i};
end

end

%% update the primal gradient

L1_v(t) = 0;
i=1;
L1_vp{i}(t) = 0;
g1{i} = zeros(size(St_im{i}));
% 
for k=1:Ps
[idx, dual_var1_, g1_, norm1_] = fetchNext(f1);
dual_var1{i,idx} = dual_var1_;
g1{i} = g1{i} + g1_;
L1_vp{i}(t) = L1_vp{i}(t) + norm1_;
end
g1{i} = sigma1(1) * g1{i};

if P>1
i =2;
L1_vp{i}(t) = 0;
g1{i} = zeros(size(St_im{i}));
for k=1:Ps
    [idx, dual_var1_, g1_, norm1_] = fetchNext(f2);
    dual_var1{i,idx} = dual_var1_;
    g1{i} = g1{i} + g1_;
    L1_vp{i}(t) = L1_vp{i}(t) + norm1_;
end
g1{i} = sigma1(1) * g1{i};

i=3;
L1_vp{i}(t) = 0;
g1{i} = zeros(size(St_im{i}));
for k=1:Ps
    [idx, dual_var1_, g1_, norm1_] = fetchNext(f3);
    dual_var1{i,idx} = dual_var1_;
    g1{i} = g1{i} + g1_;
    L1_vp{i}(t) = L1_vp{i}(t) + norm1_;
end
g1{i} = sigma1(1) * g1{i};
% end

for i=3 %:P
L1_v(t) = L1_v(t) + L1_vp{i}(t);
end



clear idx  g1_ norm1_;
%%%%%%%

%Free memory
g2_val = [];
clear ns;
%     clear prev_xsol;
St_im_old = [];
tm = toc(tm);


if method ~= 0
%**************************************************************************
% Checking the polarization constraint
%**************************************************************************

if stokes_P>1
v12 = sqrt(abs(St_im{2}).^2 + abs(St_im{3}).^2);
[c_f(t),~] = size(find(v12 - St_im{1} > 0));
c_m(t) = max(max(v12-St_im{1}));
c_mn(t) = max(max((v12-St_im{1})./St_im{1}));

% if method > 1
% 
%     % To check the constraint for epi_h1
%     [c_h1(t),~] = size(find(-St_im{1}>St_im{4}));
% 
%     % To check the constraint for epi_h2
%     [c_h2(t),~] = size(find(v12-St_im{5} > 1e-5));
% 
% 
%     if param.verbose == 2
%         if c_h1(t) > 0
%             fprintf('Epi_h1 not satisfied:%d \n',c_h1(t));
%         end
% 
%         if c_h2(t) > 0
%             fprintf('Epi_h2 not satisfied:%d \n',c_h2(t));
%         end
% 
%         if c_f(t) > 0
%             fprintf('Polarization constraint not satisfied:%d \n',c_f(t));
%         end
% 
%     end
% 
%     if (param.verbose == 2)
%         figure(102)
%         subplot 221, plot(c_h1), title('Count of unsatisfied epi h1')
%         subplot 222, plot(c_h2), title('Count of unsatisfied epi h2')
%         subplot 223, plot(c_f), title('Count of unsatisfied total flux')
%         subplot 224, plot(c_m), title('Max of norm(Q,U)-I')
%     end
% end
end
end

%% Stopping criterion and logs
% fprintf('Time for iteration %i: %3.3f\n\n\n',t, tm);

% figure(200),
% subplot 241, semilogy(L2_vv{1})
% subplot 244, semilogy(L1_v)
% subplot 245, semilogy(rel_sol_norm_change{1})
% % if stokes_P>1
%     subplot 242, semilogy(L2_vv{2})
%     subplot 243, semilogy(L2_vv{3})
%     subplot 246, semilogy(rel_sol_norm_change{2})
%     subplot 247, semilogy(rel_sol_norm_change{3})
%     
% % end

if (rem(t,200) == 0)
figure(200)
% if method ~= 0
% subplot 331, semilogy(L1_v)
% end
subplot 334, semilogy(rel_sol_norm_change{1})
subplot 335, semilogy(rel_sol_norm_change{2})
subplot 336, semilogy(rel_sol_norm_change{3})
subplot 337, plot(snr{1})
subplot 338, plot(snr{2})
subplot 339, plot(snr{3})
title(sprintf('method=%d',method))
pause(0.1)
end

im = x_th.im_true;
% 
%     figure(101)
%     subplot 231, imagesc(im{1}), title('True I'), axis image, axis off
%     subplot 232, imagesc(im{2}), title('True Q'), axis image, axis off
%     subplot 233, imagesc(im{3}), title('True U'), axis image, axis off
%     subplot 234, imagesc(St_im{1}), title('Reconstructed I'), axis image, axis off
%     subplot 235, imagesc(St_im{2}), title('Reconstructed Q'), axis image, axis off
%     subplot 236, imagesc(St_im{3}), title('Reconstructed U'), axis image, axis off

if stokes_P>1 
%  count_flux_thresh;
 flux_thresh = 0;
f_thresh(t) = flux_thresh;
end

%% Reweighting procedure

% 
% if stokes_P > 1
%     check_norm = (norm2{1,1}) <= (epsilonts{1}{1}) && (norm2{2,1}) <= (epsilonts{2}{1}) && ...
%         (norm2{3,1}) <= (epsilonts{3}{1}) && (norm2{4,1}) <= (epsilonts{4}{1}) ;
% %     if method == 3
% %         check_norm = check_norm && (f_thresh(t)<= param.pol_thresh);
% %     end
% else
%     check_norm =  (norm2{1,1}) <= (epsilonts{1}{1});
% end
% 
% 
% if (param.use_reweight_steps || param.use_reweight_eps) && ( t < param.reweight_max_reweight_itr) || (param.init_rw)
% 
% 
%     if (param.use_reweight_steps && t == reweight_steps(reweight_step_count)) || ...
%             (param.use_reweight_eps && ...
%             check_norm && ...
%             param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && ...
%             rm(t) < param.reweight_rel_obj) || ...
%             (param.init_rw && param.flag_change_method && t == param.init_stop_iter)
%         %        norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts)) && ...
% 
%         for i =1:stokes_P
%             for k = 1:Ps
%                 d_val = abs(Psit{k}(reshape(St_im{i},Ny,Nx)));
%                 d_val_max = max(d_val(:));
%                 if d_val_max == 0
%                     d_val_max =1;
%                 end
%                 weights{i,k} = param.reweight_alpha(k) ./ (param.reweight_alpha(k) + (d_val./d_val_max));
%                 weights{i,k}(d_val > d_val_max * param.reweight_abs_of_max) = 0;
%                 eta{i,k} = param.eta_o(i)*weights{i,k};
%                 param.reweight_alpha(k) = param.reweight_alpha_ff(k) .* param.reweight_alpha(k);
%             end
%         end
%         param.weights = weights;
%         fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);
% 
%         if reweight_step_count > param.total_reweights
%             param.reweight_max_reweight_itr = t+1;
%             fprintf('\n\n\n\n\n\n\n No more reweights \n\n\n\n\n');
%         end
% 
%         reweight_step_count = reweight_step_count + 1;
%         reweight_last_step_iter = t;
%         reweight_step_count_since_best_bound_search = reweight_step_count_since_best_bound_search + 1;
%         rw = rw+1;
% %         util_create_pool(9);
% 
%         
%         %%%%%%%% Added by jbirdi
%         if reweight_step_count == 2 && param.flag_change_method
%             param.init_rw = 0;
%             method = 3;
%              param.tau = 0.33;
%              
%         end
%         %%%%%%%
% 
% %         St_im_sol{rw} = St_im;
% %         dual_var1_sol{rw} = dual_var1;
% % 
% %         dual_var2_sol{rw} = dual_var2;
% % 
% %         dual_var3_sol{rw} = dual_var3;
% % 
% %         snr_v_sol{rw} = snr_v;
% % 
% %         t_sol{rw} = t;
% % 
% %         rel_sol2{rw} = rel_sol;
% % 
% %         mse_sol{rw} = mse;
% % 
% %         nrmse_sol{rw} = nrmse;
% % 
% %         f_thresh_sol{rw} = f_thresh;
% % 
% %         cf_sol{rw} = c_f;
% % 
% %         fitswrite(St_im{1},[param.path,'rwsol_I_',num2str(rw),'.fits']);
% %         fitswrite(St_im{2},[param.path,'rwsol_Q_',num2str(rw),'.fits']);
% %         fitswrite(St_im{3},[param.path,'rwsol_U_',num2str(rw),'.fits']);
% 
%     end
% end
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global stopping criteria
if rm(t) < param.rel_obj %% && f_thresh(t) <= param.pol_thresh
%     && reweight_step_count >= param.total_reweights 
%     &&  ((param.global_stop_bound && check_norm))

    %             ((param.global_stop_bound && norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts))) || ...
    %             (~param.global_stop_bound && prod(cell2mat(norm2) <= cell2mat(epsilonts))))
    %
    flag = 1;
    break;
end
end


% disp('***********************')
% % fprintf('Seed %i\n',n_test);
% fprintf('Iter %i\n',t);
% fprintf('L1 norm              = %e\n',L1_v(t))
% % fprintf('Residual             = %e\n',residual(t))
% % fprintf('Obj function value   = %e\n',obj(t))
% fprintf('Rel sol change = %e\n', rm(t))
% % fprintf('Residual = %e\n',residual(t))
% % fprintf('Flux constraint not satisfied:%d \n',c_f(t));
% fprintf('Flux constraint not satisfied after threshold:%d \n',f_thresh(t));
% % fprintf('SNR_I = %e, SNR_Q = %e, SNR_U = %e \n',snr_v{1}(t), snr_v{2}(t), snr_v{3}(t));
% 
% disp('***********************')

param_dual.dual_var1 = dual_var1;
param_dual.dual_var2 = dual_var2;
param_dual.dual_tilde2 = dual_tilde2;
param_dual.St_im = St_im; 
% param_dual.St_im{4} = St_im{4};
% param_dual.St_im{5} = St_im{5};
l1 = L1_v(end);

end

