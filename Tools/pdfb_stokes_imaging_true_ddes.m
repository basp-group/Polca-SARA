% function [St_im, L1_v, L1_vp, L2_vv, L2_vp, rel_sol_norm_change, sol_v, snr_v, no_sub_itr_v, dual_var1, dual_var3, sol_best_bound_v, sol_reweight_v] = pdfb_par_rescaled_precond_wave_par_stokes_rw(y, epsilont, epsilonts, A, At, T, pU, W, Psi, Psit, Psiw, Psitw, Gw,P, param)


%% Functions used in proximity operators

% scaling, projection on L2 norm
sc = @(z, radius) z * min(radius/norm(z(:)), 1);

% thresholding negative values
negt = @(z) max(real(z), 0);

%soft thresholding operator
soft = @(z, T) sign(z) .* max(abs(z)-T, 0);

% Positive soft thresholding operator
pos_soft = @(z, T) max(real(z)-T, 0);

% phi__= 0.5*(sqrt(5)-1); % Golden ratio
phi__ = param.phi_val;

count_eps_update_up = 0; 
count_eps_update_down = 0; 

method = param.method;
R=1;
% for i=1:P
% for q = 1:R 
%     t__{i,q}=1;
% end
% end


% % number of nodes
% R = length(y{1,1});
% pU = aW;

% T = G_;

Nx = param.Nx;
Ny = param.Ny;
Ps = param.Ps;
method = param.method;
evl = param.evl;

 def_var_real_data;
    
   % 
    init_var_real_data;
   
    
    if param.nnls_init
        for i=1:P
        St_im{i} = result_st_nnls.sol{i};
        for k = 1:param.Ps
                d_val = abs(Psit{k}(reshape(St_im{i},Ny,Nx)));
                d_val_max = max(d_val(:));
                if d_val_max == 0
                    d_val_max =1;
                end
                weights{i,k} = param.reweight_alpha(k) ./ (param.reweight_alpha(k) + (d_val./d_val_max));
                weights{i,k}(d_val > d_val_max * param.reweight_abs_of_max) = 0;
                eta{i,k} = param.eta_o(i)*weights{i,k};
                param.reweight_alpha(k) = param.reweight_alpha_ff(k) .* param.reweight_alpha(k);
            end
        end
    end
    
    
   % % Change the cell formatting
for i=1:P
     for q=1:R
         dual_var3_{i}{q} = dual_var3{i,q};
%          epsilont_{i}{q} = epsilont{i,q};
% %       
     end
end

 dual_var3 = dual_var3_;
 
 %%%%
% 
% for i=1:P
%     
% [proj{i}] = solver_find_elipse_point(y{i}, pU, A, T{i}, St_im{i}, dual_var3{i}, W{i}, epsilont{i}, param.elipse_proj_max_iter, param.elipse_proj_min_iter, param.elipse_proj_eps);
%    
% end
% % 
% A = afclean(A);
% At = afclean(At);


parfor k = 1:Ps
    Psi{k} = afclean(Psi{k});
    Psit{k} = afclean(Psit{k});
end
Psiw = [];
Psiwt = [];
util_create_pool(Ps);



% if W{1}{1} ~= ':'
%     % oversampling vectorized data length
%     No = size(W{1}{1}, 1);
% else
%     
%     No = size(T{1}{1}' * y{1}{1} , 1);
% end

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
if ~isfield(param, 'nu3')
    param.nu3 = zeros(R, 1);
    % maximum eigenvalue of operato A^T A
    for q = 1:R
        Tw = spalloc(size(T{q}, 1), No, size(T{q}, 2) * 16);
        Tw(:, W{q}) = T{q};
        fprintf('\nComputing operator norm: block %i \n', q)
        param.nu3(q) = op_norm(@(x) Tw * A(x), @(x) At(Tw' * x), [Ny, Nx], 1e-4, 200, 1);
        clear Tw;
    end
    fprintf('\n');
else
    if numel(param.nu3) == 1
        param.nu3 = ones(R, 1) * param.nu3;
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
if method == 3
    param.tau = 0.33;
end

 if ~isfield(param, 'weights')
    param.weights = cell(Ps, 1);
    for i=1:P
    for k = 1:Ps
        param.weights{i,k} = ones(size(Psit{k}(At(zeros(No, 1))), 1), 1);
    end
    end
else
    if ~iscell(param.weights)
        weights = param.weights;
        param.weights = cell(Ps, 1);
        for i=1:P
        for k = 1:Ps
            param.weights{i,k} = weights;
        end
    end
    end
 end

%  if r == 1 || rw == 0
if ~param.nnls_init
     for i=1:P
         for k=1:Ps
             param.weights{i,k} = ones(size(Psit{k}(At(zeros(No, 1))), 1), 1);
             eta{i,k} = param.eta_o(i)*param.weights{i,k};
         end
     end
end
%  end

 %% set up log variables
L1_v = zeros(param.max_iter, 1);
% L1_vp = zeros(param.max_iter, Ps);
% L2_v = zeros(param.max_iter, 1);
% L2_vp = zeros(param.max_iter, R);
% no_sub_itr_v = cell(param.max_iter, 1);

delta_v = zeros(param.max_iter, 1);

sol_steps = param.sol_steps;
sol_step_count = 1;

reweight_steps = param.reweight_steps;
reweight_step_count = 1;
reweight_last_step_iter = 1;
reweight_step_count_since_best_bound_search = 0;

% best_bound_steps = param.best_bound_steps;
% best_bound_step_count = 1;

if (param.use_best_bound_steps ||param.use_best_bound_eps)
   sol_best_bound_v = zeros(0, Ny, Nx);
else
   sol_best_bound_v =[];
end

sol_v =[];% zeros(length(sol_steps)-1, Ny, Nx);

sol_reweight_v = [];% zeros(0, Ny, Nx);

snr_v = [];%zeros(param.max_iter, 1);

%% store useful variables
% step size for the dual variables
sigma1 = param.sigma1;
sigma2 = param.sigma2;
sigma3 = param.sigma3;

% step size primal 
tau = param.tau;

% relaxation parameters
lambda0 = param.lambda0;
lambda1 = param.lambda1;
lambda2 = param.lambda2;


% weights
weights = param.weights;
param.weights = [];


reweight_steps = param.reweight_steps;
reweight_step_count = 1;
reweight_last_step_iter = 1;
reweight_step_count_since_best_bound_search = 0;
rw = 0;


 %%
%****************************************************************************************%
%****************************************************************************************%

% Algorithm

g = zeros(size(St_im{1})); %gradient part 
use_proj_elipse_fb =  param.use_proj_elipse_fb;
elipse_proj_max_iter = param.elipse_proj_max_iter; 
elipse_proj_min_iter =param.elipse_proj_min_iter;
elipse_proj_eps = param.elipse_proj_eps;

% Tt = cell(P, R);
% for i=1:P
% for q = 1:R
%     Tt{i,q} = T{i}{q}';
% end

% Initialize variables 
for i=1:stokes_P
St_im_old{i} = St_im{i};
g1{i} = zeros(size(St_im{i}));
g2{i} = zeros(size(St_im{i}));
g3{i} = zeros(size(St_im{i}));
prev_xsol{i} = 2*St_im{i}-St_im_old{i};
end


%%

for t = 1:param.max_iter
tm = tic;

%% primal update

if mod(t,1000) ==0
    % Restored images
    fitswrite(St_im{1},[param.path,'solI_',num2str(t),'.fits']);
    fitswrite(St_im{2},[param.path,'solQ_',num2str(t),'.fits']);
    fitswrite(St_im{3},[param.path,'solU_',num2str(t),'.fits']);
    % Residual images
    fitswrite(real(res_im{1}),[param.path,'resI_',num2str(t),'.fits']);
    fitswrite(real(res_im{2}),[param.path,'resQ_',num2str(t),'.fits']);
    fitswrite(real(res_im{3}),[param.path,'resU_',num2str(t),'.fits']);
end


for i=1:P

%% Primal update
St_im_old{i} = St_im{i};
if method == 1
    St_im{i} = real((St_im{i} - tau*(reshape(g1{i},Ny,Nx) + reshape(g3{i},Ny,Nx))));
else
    St_im{i} = real((St_im{i} - tau*(reshape(g1{i},Ny,Nx) + reshape(g2{i},Ny,Nx)+reshape(g3{i},Ny,Nx))));
end
if i == 1
    St_im{i} = negt(St_im{i});
end
prev_xsol{i} = 2*St_im{i}-St_im_old{i};


%% Relative change in the estimates
norm_prevsol(i) = norm(St_im_old{i}(:));
rm(t) = 0;
if (norm_prevsol(i) == 0)
    rel_sol_norm_change{i}(t) = 1;
else
    rel_sol_norm_change{i}(t) = norm(St_im{i}(:) - St_im_old{i}(:))/norm(St_im{i}(:));
    rm(t) = max(rel_sol_norm_change{i}(t),rm(t));
end

end

%% Dual variables update

%% L1 prox update
%%%%%#### Parallelize for the sparsity bases

i=1;

L1_vp{i}(t) = 0;
g1{i} = zeros(size(St_im{i}));
for k = 1:Ps
f1(k) = parfeval(@run_par_waverec2,3,dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{i}, eta{i,k});
end

if P>1
i=2;
for k = 1:Ps
    f2(k) = parfeval(@run_par_waverec2,3,dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{i}, eta{i,k});
end

i=3;
for k = 1:Ps
    f3(k) = parfeval(@run_par_waverec2,3,dual_var1{i,k},Psit{k}, Psi{k}, prev_xsol{i}, sigma1(k), St_im{i}, eta{i,k});
end
end


%% L2 ball projection update

residual(t) = 0;

if blocks
for i=1:P

% non gridded measurements of current solution
ns{i} = A(prev_xsol{i});
g2_val = zeros(No,1);

% select parts to be sent to nodes
for q = 1:R

    r2{i,q} = T{i}{q} * ns{i}(W{i}{q});

    if use_proj_elipse_fb
        [proj{i}{q}, no_sub_itr{i}{q}] = solver_proj_elipse_fb(1 ./ pU{1}  .* (dual_var3{i}{q}), r2{i,q}, y{i}{q}, pU{1}, epsilont{i}{q}, proj{i}{q}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
        dual_var3{i}{q} = dual_var3{i}{q} + pU{1} .* (r2{i,q} - proj{i}{q});
    end
    %         if ~param.use_proj_elipse_fb
    %             [proj{i,q}, no_sub_itr{i,q}] = solver_proj_elipse(1 ./ sqrt(pU{i,q}) .* (dual_var3{i}{q}./sigma3(i) + pU{i,q} .* r2{i,q}), 1 ./ sqrt(pU{i,q}), epsilont{i,q}, proj{i,q}, y{i}{q}, elipse_proj_max_iter, elipse_proj_min_iter, elipse_proj_eps);
    %             vy2{i,q} = dual_var3{i}{q} + pU{i,q} .* r2{i,q} - sqrt(pU{i,q}) .* proj{i,q};
    %         end

    u2{i,q} = Tt{i,q} * dual_var3{i}{q};
    g2_val(W{i}{q}) = g2_val(W{i}{q}) + u2{i,q};

    % norm of residual
    res{i,q} = r2{i,q} - y{i}{q};
    norm2{i,q} = norm(res{i,q});
    residual(t) = residual(t) + norm2{i,q};
end

L2_vp{i}(t) = norm2{i,q};

g3{i} = sigma3* (At(g2_val));
no_sub_itr_v{i}(t) = no_sub_itr{i};


%%%% Residual images

if mod(t,999) == 0
    res_im{i} = At(Gw{i}' * (y{i}{1} - Gw{i} * A(St_im{i})));
end
end
else
    b = conv_stokes_to_bright([prev_xsol{1}(:),prev_xsol{2}(:),prev_xsol{3}(:)],L,size(prev_xsol{1}));
    for i = 1:4
    
     r2{i,q} = A{i,1}(b{i});
     
     %%%%%%%%
     if sep == 1
            dual_tilde3{i} = dual_var3{i} + r2{i,q};
            b1{i} = sc((dual_tilde3{i})-y(:,i),epsilon{i})+y(:,i);
            dual_var3{i} = dual_tilde3{i} - b1{i};
     else 
         
       
            dual_tilde3(:,i) = dual_var3(:,j) + sigma3*phi_s{j};
            
        
            
 
        dt3 = dual_tilde3(:);
        b1 = sc((dt3./sigma3) - y_conc(:), epsilon{n_test}) + y_conc(:);
        bb = vec2mat(b1',5631);
        bb = bb';
        dual_var3 = dual_tilde3 - sigma3*bb;
     end
    
     
     % norm of residual
    res{i,q} = r2{i,q} - y(:,i);
    norm2{i,q} = norm(res{i,q});
    residual(t) = residual(t) + norm2{i,q};

L2_vp{i}(t) = norm2{i,q};
g3_{i} = (At{i,1}(dual_var3{i}));
    end


% g3{i} = sigma3* (At{i,1}(dual_var3{i}));

g3 = sigma3*mat2cell([g3_{1:4}]*Lt,size(g3_{1},1), [1 1 1])';

% no_sub_itr_v{i}(t) = no_sub_itr{i};
end

          %%%%%%%%
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
for i=1:P+2
    dual_tilde2{i} = dual_var2{i} + prev_xsol{i};
end

% par
vec = [dual_tilde2{1}(:),dual_tilde2{4}(:)];
epi_h{1} = vec - Proj_epih1(1,vec);

vec = [dual_tilde2{2}(:),dual_tilde2{3}(:),dual_tilde2{5}(:)];
epi_h{2} = vec - Proj_epih2(1,vec);


% Update variables
dual_var2{1} = reshape(epi_h{1}(:,1),Ny,Nx);
dual_var2{4} = reshape(epi_h{1}(:,2),Ny,Nx);

dual_var2{2} = reshape(epi_h{2}(:,1),Ny,Nx);
dual_var2{3} = reshape(epi_h{2}(:,2),Ny,Nx);
dual_var2{5} = reshape(epi_h{2}(:,3),Ny,Nx);

for i=1:P+2
    g2{i} = sigma2*dual_var2{i};
end

end


%Free memory
g2_val = [];
clear ns;
%     clear prev_xsol;
St_im_old = [];


%% To estimate the epsilon bounds
% 
% 
%     %%##### AR ADAPTIVE bound update on each block
% if( param.use_adapt_bound_eps ==1)  %%##### AR
% 
%     for i = 1:P
%         for q = 1:R
%                         if (norm2{i,q}<epsilonts{i}{q}) && (norm2{i,q}> (1-param.adapt_bound_tol)*epsilont{i}{q})
%                  t__{i,q} = param.max_iter;
%                         elseif   ((t>t__{i,q}+param.adapt_bound_steps) || (t__{i,q} == param.max_iter )) && ...
%                     (norm2{i,q}< (1-param.adapt_bound_tol)*epsilont{i}{q}) &&  ...
%                     (rel_sol_norm_change{i}(t) <  param.reweight_rel_obj) 
%                 t__{i,q} = t;
% 
%                 epsilont{i}{q} = norm2{i,q} + (-norm2{i,q} + epsilont{i}{q})*(1-phi__);
%                 epsilonts{i}{q}= 1.001*epsilont{i}{q};
%                 %                      epsilon = norm(cell2mat(epsilont{i,:}));
%                 %                      epsilons= 1.001*epsilon;
%                 count_eps_update_down = count_eps_update_down +1;
%                 disp (['#####     Updated  epsilon DOWN: ', num2str(epsilont{i}{q})])
%                 disp ''
%             end
% 
%             if (norm2{i,q}>epsilonts{i}{q}) 
%                 if  (( (t>t__{i,q}+param.adapt_bound_steps) && (rel_sol_norm_change{i}(t) < param.adapt_bound_rel_obj)) || ...
%                         ((t__{i,q} == param.max_iter ) && (rel_sol_norm_change{i}(t) <  param.reweight_rel_obj))    ) || ...
%                        ( ( t ==param.adapt_bound_start)) 
% %                    || (t > t__{i,q}+param.adapt_steps))
%                     t__{i,q} = t;
%                     %                         disp ''
% %                     disp (['#####    Updated bound-', num2str(count_eps_update_up)])
%                     epsilont{i}{q} =epsilont{i}{q} + (norm2{i,q} - epsilont{i}{q})*phi__;% (norm2{q} + epsilont{q,p})/2;
% %%%%% Addedd by jbirdi
% % epsilont{i}{q} = epsilont{i}{q}*(1+param.eps_update_tol);
% %%%%%%%
%                     epsilonts{i}{q}= 1.001*epsilont{i}{q};
%                     %                         epsilon = norm(cell2mat(epsilont{i,:}));
%                     %                         epsilons= 1.001*epsilon;
%                     count_eps_update_up = count_eps_update_up +1;
%                     disp(['#####     Updated  epsilon UP: ', num2str(epsilont{i}{q})])
%                     disp ''
%                 end
%             
%             end
%         end
%     end
% end

    %##### AR ADAPTIVE bound update on each block

% 
% if param.use_adapt_bound_eps == 1
% for i=1:P
%     for q=1:R
%         if norm2{i,q} < (1-param.adapt_bound_tol)*epsilont{i}{q}
%             if t > t__{i,q}+param.adapt_bound_steps && rel_sol_norm_change{i}(t) < param.adapt_bound_rel_obj
%                 t__{i,q} = t;
%                 epsilont{i}{q} = norm2{i,q} + (-norm2{i,q} + epsilont{i}{q})*(1-phi__);
%                 count_eps_update_down = count_eps_update_down +1;
%                 disp (['#####     Updated  epsilon DOWN: ', num2str(epsilont{i}{q})])
%                 disp ''
%             end
%         end
% 
%         if norm2{i,q} > (1+param.adapt_bound_tol)*epsilont{i}{q}
%             if t > t__{i,q}+param.adapt_bound_steps && rel_sol_norm_change{i}(t) < param.adapt_bound_rel_obj
%                 t__{i,q} = t;
%                 epsilont{i}{q} =epsilont{i}{q} + (norm2{i,q} - epsilont{i}{q})*phi__;% (norm2{q} + epsilont{q,p})/2;
%                 count_eps_update_up = count_eps_update_up +1;
%                 disp(['#####     Updated  epsilon UP: ', num2str(epsilont{i}{q})])
%                 disp ''
%             end
%         end
%     end
% end
% end



%% update the primal gradient

L1_v(t) = 0;
i=1;
L1_vp{i}(t) = 0;
g1{i} = zeros(size(St_im{i}));

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
%%%%%%%

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

end

for i=1:P
L1_v(t) = L1_v(t) + L1_vp{i}(t);
L2_vv{i}(t) = norm((norm2{i}));
end



clear idx  g1_ norm1_;
tm = toc(tm);



%**************************************************************************
% Checking the polarization constraint
%**************************************************************************

if P>1
v12 = sqrt(abs(St_im{2}).^2 + abs(St_im{3}).^2);
[c_f(t),~] = size(find(v12 - St_im{1} > 0));
c_m(t) = max(max(v12-St_im{1}));
c_mn(t) = max(max((v12-St_im{1})./St_im{1}));

if method > 1

    % To check the constraint for epi_h1
    [c_h1(t),~] = size(find(-St_im{1}>St_im{4}));

    % To check the constraint for epi_h2
    [c_h2(t),~] = size(find(v12-St_im{5} > 1e-5));


    if param.verbose == 2
        if c_h1(t) > 0
            fprintf('Epi_h1 not satisfied:%d \n',c_h1(t));
        end

        if c_h2(t) > 0
            fprintf('Epi_h2 not satisfied:%d \n',c_h2(t));
        end

        if c_f(t) > 0
            fprintf('Flux constraint not satisfied:%d \n',c_f(t));
        end

    end

    if (param.verbose == 2)
        figure(102)
        subplot 221, plot(c_h1), title('Count of unsatisfied epi h1')
        subplot 222, plot(c_h2), title('Count of unsatisfied epi h2')
        subplot 223, plot(c_f), title('Count of unsatisfied total flux')
        subplot 224, plot(c_m), title('Max of norm(Q,U)-I')
    end
end
end

%% Stopping criterion and logs
fprintf('Time for iteration %i: %3.3f\n\n\n',t, tm);

figure(200),
subplot 241, semilogy(L2_vv{1})
subplot 244, semilogy(L1_v)
subplot 245, semilogy(rel_sol_norm_change{1})
if P>1
    subplot 242, semilogy(L2_vv{2})
    subplot 243, semilogy(L2_vv{3})
    subplot 246, semilogy(rel_sol_norm_change{2})
    subplot 247, semilogy(rel_sol_norm_change{3})
end


%
%     for p = 1:P
%     ee(p,t) = epsilont{q,p};
%     end
%
% %     if rem(t,50) == 0
% %         figure(600),
% %         subplot 131, semilogy(ee(1,:))
% %         subplot 132, semilogy(ee(2,:))
% %         subplot 133, semilogy(ee(3,:))
% %         pause(0.1)
% %     end
%
% log
if (param.verbose >=4)
    fprintf('Iter %i\n',t);
    fprintf('L1 norm              = %e\n',L1_v(t))
    fprintf('Residual             = %e\n',residual(t))
    %         fprintf('Obj function value   = %e\n',obj(t))
    %         fprintf('Rel sol norm change  = %e\n',rel_sol)
end

if (param.verbose == 4 && rem(t,10) == 0)
    figure(100)
    subplot 331, semilogy(L1_v), title('Obj function')
    subplot 334, plot(snr_v{1}), title('SNR- I')
    subplot 335, plot(snr_v{2}), title('SNR- Q')
    subplot 336, plot(snr_v{3}), title('SNR- U')
    subplot 337, semilogy(rel_sol{1}), title('rel norm- I')
    subplot 338, semilogy(rel_sol{2}), title('rel norm- Q')
    subplot 339, semilogy(rel_sol{3}), title('rel norm- U')

    %         figure(200)
    %         subplot 121, semilogy(L1_v), title('l1 norm')
    %         subplot 122, semilogy(residual), title('residual')
    %
    %
    %
    figure(101)
    subplot 231, imagesc(im{1}), title('True I'), axis image, axis off
    subplot 232, imagesc(im{2}), title('True Q'), axis image, axis off
    subplot 233, imagesc(im{3}), title('True U'), axis image, axis off
    subplot 234, imagesc(St_im{1}), title('Reconstructed I'), axis image, axis off
    subplot 235, imagesc(St_im{2}), title('Reconstructed Q'), axis image, axis off
    subplot 236, imagesc(St_im{3}), title('Reconstructed U'), axis image, axis off

    figure(501), imagesc(mask2)

    %         figure(102)
    %         subplot 211, plot(c_f), title('Count of unsatisfied total flux')
    %         subplot 212, plot(min_res)

    pause(0.1)
end

if P>1
count_flux_thresh;
f_thresh(t) = flux_thresh;
end

%% Reweighting procedure

if P > 1
    check_norm = (norm2{1,1}) <= (epsilonts{1}{1}) && (norm2{2,1}) <= (epsilonts{2}{1}) && ...
        (norm2{3,1}) <= (epsilonts{3}{1});
    if method == 3
        check_norm = check_norm && (f_thresh(t)<= param.pol_thresh);
    end
else
    check_norm =  (norm2{1,1}) <= (epsilonts{1}{1});
end

if (param.use_reweight_steps || param.use_reweight_eps) && ( t < param.reweight_max_reweight_itr) || (param.init_rw)


    if (param.use_reweight_steps && t == reweight_steps(reweight_step_count)) || ...
            (param.use_reweight_eps && ...
            check_norm && ...
            param.reweight_min_steps_rel_obj < t - reweight_last_step_iter && ...
            rm(t) < param.reweight_rel_obj) || ...
            (param.init_rw && param.flag_change_method && t == param.init_stop_iter)
        %        norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts)) && ...

        for i =1:P
            for k = 1:Ps
                d_val = abs(Psit{k}(reshape(St_im{i},Ny,Nx)));
                d_val_max = max(d_val(:));
                if d_val_max == 0
                    d_val_max =1;
                end
                weights{i,k} = param.reweight_alpha(k) ./ (param.reweight_alpha(k) + (d_val./d_val_max));
                weights{i,k}(d_val > d_val_max * param.reweight_abs_of_max) = 0;
                eta{i,k} = param.eta_o(i)*weights{i,k};
                param.reweight_alpha(k) = param.reweight_alpha_ff(k) .* param.reweight_alpha(k);
            end
        end
        param.weights = weights;
        fprintf('\n\n\n\n\n\n\n Performed reweight no %d \n\n\n\n\n', reweight_step_count);

        if reweight_step_count > param.total_reweights
            param.reweight_max_reweight_itr = t+1;
            fprintf('\n\n\n\n\n\n\n No more reweights \n\n\n\n\n');
        end

        reweight_step_count = reweight_step_count + 1;
        reweight_last_step_iter = t;
        reweight_step_count_since_best_bound_search = reweight_step_count_since_best_bound_search + 1;
        rw = rw+1;
        util_create_pool(9);

        
        %%%%%%%% Added by jbirdi
        if reweight_step_count == 2 && param.flag_change_method
            param.init_rw = 0;
            method = 3;
             param.tau = 0.33;
             
        end
        %%%%%%%%
% 
%         St_im_sol{rw} = St_im;
%         dual_var1_sol{rw} = dual_var1;
% 
%         dual_var2_sol{rw} = dual_var2;
% 
%         dual_var3_sol{rw} = dual_var3;
% 
%         snr_v_sol{rw} = snr_v;
% 
%         t_sol{rw} = t;
% 
%         rel_sol2{rw} = rel_sol;
% 
%         mse_sol{rw} = mse;
% 
%         nrmse_sol{rw} = nrmse;
% 
%         f_thresh_sol{rw} = f_thresh;
% 
%         cf_sol{rw} = c_f;

        fitswrite(St_im{1},[param.path,'rwsol_I_',num2str(rw),'.fits']);
        fitswrite(St_im{2},[param.path,'rwsol_Q_',num2str(rw),'.fits']);
        fitswrite(St_im{3},[param.path,'rwsol_U_',num2str(rw),'.fits']);

    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global stopping criteria
if rm(t) < param.rel_obj  && reweight_step_count >= param.total_reweights && ...
        ((param.global_stop_bound && check_norm))

    %             ((param.global_stop_bound && norm(cell2mat(norm2)) <= norm(cell2mat(epsilonts))) || ...
    %             (~param.global_stop_bound && prod(cell2mat(norm2) <= cell2mat(epsilonts))))
    %
    flag = 1;
    break;
end
end


disp('***********************')
% fprintf('Seed %i\n',n_test);
fprintf('Iter %i\n',t);
fprintf('L1 norm              = %e\n',L1_v(t))
fprintf('Residual             = %e\n',residual(t))
% fprintf('Obj function value   = %e\n',obj(t))
fprintf('Rel sol change = %e\n', rm)
fprintf('Residual = %e\n',residual(t))
fprintf('Flux constraint not satisfied:%d \n',c_f(t));
fprintf('Flux constraint not satisfied after threshold:%d \n',f_thresh(t));
% fprintf('SNR_I = %e, SNR_Q = %e, SNR_U = %e \n',snr_v{1}(t), snr_v{2}(t), snr_v{3}(t));

disp('***********************')



