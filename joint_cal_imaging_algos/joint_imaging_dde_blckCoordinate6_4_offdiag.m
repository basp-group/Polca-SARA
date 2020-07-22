function [U1, U2, epsilon, objective, time_tot, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate6_4_offdiag(y, Y, x0, epsilon, x_th, u, v, Omega, na, T, J, P, F, K, N, S, U1, U2, scale, V, param_im, param_dde, param_algo, ...
     W, name, mu,param_pd,L,Lt,B_tmp_test,Bt_tmp_test, eta1,U_th)
% Joint imaging and calibration procedure in absence of temporal
% regularization.
%-------------------------------------------------------------------------%
%%
% Input:
% < y       : data vector 
% < Y       : data with redundancy
%             (Nan on the diag., y_ab = conj(y_ab) for a > b)
% < x0      : reference image
% < epsilon : initial perturbation (generally zeros(size(x0)))
% < x_th    : ground truth image
% < u       : u component [M, T]
% < v       : v component [M, T]
% < Omega   : u-v components with redundancy (similar to Y) [na, na, 2, T],
%             Omega(:, :, 1, t) -> v(t)
%             Omega(:, :, 2, t) -> u(t)
% < na      : number of antennas
% < T       : number of snapshots
% < J       : size of the gridding kernels (square kernels [J, J])
% < P       : size of the DDE support in the temporal Fourier space
% < F       : size of the temporal Fourier space (F = T)
% < K       : size of the spatial Fourier space [2, 1]
% < N       : image dimension [2, 1]
% < S       : spatial dimension of the DDE kernels (square kernels [S, S])
% < U1/U2   : initial DDE kernels (in the temporal Fourier domain) [S2, na, P]
% < scale   : scaling coefficients from the NUFFT (size(scale) = N)
% < V       : spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na, na, T]
% < param_im  : image-related parameters ()
% < param_dde : DDE-related parameters (contraints definition)
%               - .theta_maxR (.theta_maxI): maximum value for the real
%                                             (imag.) part
%               - .theta_minR (.theta_minI): minimum value for the real
%                                             (imag.) part
%               - .max_it :maximum number of iteration of the malor loop
%               - .JUtot : number of iteration for the DDEs
%               - .JU1 : number of iteration to estimate U1
%               - .JU2 : number of iteration to estimate U2
%               - .nu : hyperparameter controlling the Eucl. distance
%                        between U1 and U2
%               - .tol_norm : stopping criterion of the DDE major loop
% < param_algo : algorithm parameters
%                - .min_x / .max_x : bounds on the perturbation espilon
% < W          : spatial NUFFT gridding kernels [J2, T]
% < name       : name of the temporary backup file
% < mu         : additional hyperparameter (see eq (13) in Audrey's SPIE paper)

% Output:

% > U1/U2       : estimated DDEs [S2, na, P]   
% < epsilon     : estimated pertubration
% < objective   : objective function
% < time_tot    : computation time for each iteration of the outermost loop
% < error_dirty : SNR of the dirty image
% < snr_dirty   : SNR of the dirty image
% < snr_x       : image reconstruction quality (in dB)
%-------------------------------------------------------------------------%
%% 

%% Initialization
time_tot = zeros(param_dde.max_it+1, 1);
objective = zeros(param_dde.max_it+1, 1);
error_dirty = zeros(param_dde.max_it+1, 1);
snr_dirty = zeros(param_dde.max_it+1, 1);
snr_x = zeros(param_dde.max_it+1, 1);
M = param_pd.M;

param_l1.min_x = param_im.min_x;
param_l1.max_x = param_im.max_x;
param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
% param_l1.weights = 1;
param_l1.eta_o = eta1; %[0, eta1,0];
% param_l1.eta_o = [1e+4,1e+3,1e+3]; % [8e+8,2e+8,2e+8]; %[1e-4, 1e-2, 1e-2];
param_l1.approx = 1; %1; % If box bounds on the image to be considered
util_create_pool(param_pd.Ps);
lambda_scale = 1; %1/prod(K);
param_algo.eta = param_l1.eta_o*lambda_scale;
%Ax = 7.660669e+06 *2;
for i =1:3
    param_l1.mask_app{i} = (abs(x0{i}) > 0);
    param_l1.xA{i} = x0{i};
    if param_pd.dual_fb == 0
        for k=1:param_pd.Ps
            Psitw{k} = param_l1.Psit{k};
            Psiw{k} = param_l1.Psi{k};
            x_th.psit_xo{i,k} = Psitw{k}(x_th.xo{i});
            param_l1.weights{i,k} = ones(size(param_l1.Psit{k}(x_th.xo{i})));
            eta{i,k} = param_l1.eta_o(i)*lambda_scale*param_l1.weights{i,k}; %*1.9/Ax;

        end
    else
                param_l1.weights{i} = ones(size(param_l1.Psit(x_th.xo{i})));
                eta = param_algo.eta; 

    end
end



proxEpsilon = @(x,Ax) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1);
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));
if param_pd.dual_fb
    regul_x = @(x,x0,k,i) param_algo.eta(i) * sum( abs(param_l1.Psit(x0+x)) ) ;
else
    regul_x = @(x,x0,k,i) param_algo.eta(i) * sum( abs(param_l1.Psit{k}(x0+x)) ) ;
end

theta_maxoR = param_dde.theta_maxR;
theta_minoR = param_dde.theta_minR;
theta_maxoI = param_dde.theta_maxI;
theta_minoI = param_dde.theta_minI;
JU2o = param_dde.JU2;
JU1o = param_dde.JU1;
nuo = param_dde.nu;
 Jeps = param_algo.Jeps;
tol_x = param_algo.tol_x;


% Create Gt operator
Gt = []; % unsued parameter (left for legacy purpose, to be cleansed later...)
scale_t = [];
S2 = S^2;

% Initial value of the objective function
G = cell(T, 1);
[~,D1] = computeD_reg(U1, T, Gt, scale_t, T);
% D2 = computeD_reg(U2, T, Gt, scale_t, T);
D2 = D1;

% 
% 
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
    B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, scale);
    Bt_tmp{i} = @(x) real(so_fft2_adj((G0_{i})'*x, N, K, scale));
end

% Phi = zeros(size(U1));
% Phi(floor(S2/2)+1, :, :,floor(P/2)+1) = 1;

for i =1:3
    x{i} = x0{i} + epsilon{i};
end

 b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));

%  b = x_th.bright_im;

% 
% data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
% err1 = nuo*sum(abs(U1(:) - U2(:)).^2)/2;
% err2 = mu*(sum(abs(U1(:) - Phi(:)).^2) + sum(abs(U2(:) - Phi(:)).^2))/2; 
% objective(1) = data_fid(1) + err1 + err2 + regul_x(epsilon);
% error_dirty(1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
% snr_dirty(1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
% snr_x(1) = SNR(x0 + epsilon, x_th);

%% Global algorithm
fprintf('================================================================================\n');
fprintf('Joint calibration and imaging (temporal regularization) [%s]\n', name);
fprintf('================================================================================\n');
tic
for nIter = 1:param_dde.max_it
    %%% Update DDEs
    for i=1:4
        x_hat{i} = reshape(fftshift(reshape(so_fft2(b{i}, K, scale), K)), [prod(K), 1]); % FT of brightness matrix
    end
    
    % Computation of the portions of x_hat to be transmitted to the workers
    id = true(na, 1);
    X1 = cell(na, 4);
    X2 = cell(na, 4);
    for a = 1:na % relatively time consuming -> apparently no viable alternative for the later data distribution, no parfor allowed here (memory considerations)
        id(a) = false;
        idx1 = [];
        idx2 = [];
        for t = 1:T
            % for U1
            id_a = find(~isnan(squeeze(Y(a, id, t, 1))));
            om1 = squeeze(Omega(a, id, 1, t)).';
            om1 = om1(id_a);
            om2 = squeeze(Omega(a, id, 2, t)).';
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); 
            idx1 = unique([idx1; idt]); % [S2^2*J2*(na-1), na]
            % for U2
            id_a = find(~isnan(squeeze(Y(id, a, t))));
            om1 = squeeze(Omega(id, a, 1, t)); % check dimensions
            om1 = om1(id_a);
            om2 = squeeze(Omega(id, a, 2, t));
            om2 = om2(id_a);
            idt = unique(indices4Xhat(S, J, K, [om1, om2])); 
            idx2 = unique([idx2; idt]);
        end
        id(a) = true;
        for i = 1:4
            X1{a,i} = sparse(idx1, ones(size(idx1)), x_hat{i}(idx1), prod(K), 1); 
            X2{a,i} = sparse(idx2, ones(size(idx2)), x_hat{i}(idx2), prod(K), 1);
        end
    end
    clear om1 om2 x_hat
    
    for nIter_dde = 1:param_dde.JUtot
        fprintf('--DDE iterations: %4i / %4i \n', nIter_dde, param_dde.JUtot)
        
        % Update U1 & U2
        U_old = U1;
%         par
        tm = tic;
%         profile on
        %par
      for a =1:na % see parfeval for the update of U1 and U2: compute indices, then send appropriate portion of x_hat
            % update U1a
            id = true(na, 1);
            id(a) = false;
            [Ya] = create_Ya1(Y, id,a);
            Va = squeeze(V(:, a, :, :));
            Omega_a = squeeze(Omega(a,:,:,:)); % [na, 1, 2, T] om1 = squeeze(Omega(a,id,1,t)).';
            % build H1_a [na-1, S2, T], n = na-1
            Ha = zeros(8, na-1, S2, T);
            for t = 1:T % most time consuming operation
                id_a = find(~isnan(Ya{a}{1}(:, t))); % find Y(a, id, t)
                
                
                om1 = squeeze(Omega_a(id,1,t));
                om1 = om1(id_a);
                om2 = squeeze(Omega_a(id,2,t));
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); 
                if ~isempty(id_a)
                    D_at = D2{t}(:,:,id); 
                    
                     [Ha(1,id_a,:,t),Ha(2,id_a,:,t),Ha(3,id_a,:,t),Ha(4,id_a,:,t),Ha(5,id_a,:,t),Ha(6,id_a,:,t),Ha(7,id_a,:,t),Ha(8,id_a,:,t)]  = createH1a_32_reg(full(X1{a,1}(idt)), full(X1{a,2}(idt)), full(X1{a,3}(idt)), full(X1{a,4}(idt)), ...
                         Va(:, id, t), D_at, J, S, id_a); % [n, S2]
                end
                
            end

            
            [U1(:, a, :,:)] = updateUa1221_reg_offdiag(Ya{a}, reshape(U1(:, a,:, :), [S2, P, 4]), ...
                reshape(U2(:, a, :,:), [S2,P,4]), Ha, Gt, scale_t, F, T, JU1o, nuo, mu, ...
                lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            
        end

        % update D1
        [~,D1] = computeD_reg(U1, T, Gt, scale_t, T);
        stop_dde = norm(U1(:) - U_old(:))/norm(U_old(:));
                stop_dde1 = norm(U1(:) - U_th(:))/norm(U_th(:));

                fprintf('U1 estimation relative variation:%3.3f\n', stop_dde1);

%          tm = toc(tm);
%     
%     fprintf('Time elapsed:%3.3f\n',tm);
%         profile viewer
        
        % update U2a
%         U1 = U2;
        U_old = U2;
%         par
%        tm = tic;
       for a = 1:na   
            id = true(na, 1);
            id(a) = false; 
           [Ya] = create_Ya2(Y, id,a);
             
            Va = squeeze(V(:, :, a, :));
            Omega_a = squeeze(Omega(:,a,:,:)); % [na, 1, 2, T]
            % build H2_a, [na-1, S2, T], n = na-1
            Ha = zeros(8, na-1, S2, T);
            for t = 1:T
                id_a = find(~isnan(Ya{a}{1}(:, t)));
                om1 = squeeze(Omega_a(id,1,t));
                om1 = om1(id_a);
                om2 = squeeze(Omega_a(id,2,t));
                om2 = om2(id_a);
                idt = indices4Xhat(S, J, K, [om1, om2]); % [S2^2*J2*(na-1), 1] % /!\ nonzeros systematically returns a column vector
                if ~isempty(id_a)
                    D_at = D1{t}(:, :, id); % [S2, na-1, 1]
                    [Ha(1,id_a,:,t),Ha(2,id_a,:,t),Ha(3,id_a,:,t),Ha(4,id_a,:,t),Ha(5,id_a,:,t),Ha(6,id_a,:,t),Ha(7,id_a,:,t),Ha(8,id_a,:,t)] = createH2a_32_reg(full(X2{a,1}(idt)), full(X2{a,2}(idt)), full(X2{a,3}(idt)), full(X2{a,4}(idt)), ...
                        squeeze(Va(:, id, t)), D_at, J, S, id_a); % -> revoir format d'entr�e + Va
                end
            end
            U2(:, a,:, :) = updateUa2221_reg_offdiag(Ya{a}, reshape(U1(:, a, :,:), [S2, P,4]), ...
                          reshape(U2(:, a, :,:), [S2, P,4]), Ha, Gt, scale_t, F, T, JU2o, nuo, mu, ...
                          lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
        end
        % update D2
        [~,D2] = computeD_reg(U2, T, Gt, scale_t, T);
       
        % DDE stopping criterion
%         stop_dde2 = norm(U2(:) - U_old(:))/norm(U_old(:));
        stop_dde2 = norm(U2(:) - U_th(:))/norm(U_th(:));

        stop_dde = max(stop_dde, norm(U2(:) - U_old(:))/norm(U_old(:)));
        stop_crit = max(stop_dde);
        fprintf('U2 estimation relative variation:%3.3f\n', stop_dde2);
        if stop_crit(end) < param_dde.tol_norm
            fprintf('--DDE iterations: stopping criterion satisfied');
            break
        end      
        
         
    tm = toc(tm);
    
    fprintf('Time elapsed:%3.3f\n',tm);
    
    end 
    fprintf('--------------------------------------------------------------------------------\n');
    clear X1 X2
   
    %% Update Image
    
    %Create convolution matrix G
    tm = tic;
   
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
     tm = toc(tm);
    
  
    
    for i = 1:4
        B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, scale);
        Bt_tmp{i} = @(x) real(so_fft2_adj((G0_{i})'*x, N, K, scale));
    end
%    
    for i =1:3
        x{i} = x0{i} + epsilon{i};
        
    end
   
   
    b_eps = conv_stokes_to_bright([epsilon{1}(:),epsilon{2}(:),epsilon{3}(:)],L,size(epsilon{1}));
%     [B_eps,~] = Phi(B_tmp, b_eps, M);
    
    b_xo = conv_stokes_to_bright([x0{1}(:),x0{2}(:),x0{3}(:)],L,size(x0{1}));
    [~,Y_xo] = Phi(B_tmp, b_xo, M);
    
    %%%%%%%%%%
  
    
    %%%%%%
    if nIter == 1
        param_dual.dual_vv = 0;
    end
    
    [epsilon, param_dual] = updateEpsilon2_stokes(y, Y_xo, epsilon, Jeps, tol_x, nIter, proxEpsilon, B_tmp, Bt_tmp, lambda_scale,L, Lt, M, eta, param_l1,param_pd,x_th, param_dual);
    reg = 0;
    
    for i =1:3
        for k =1:param_pd.Ps
        x{i} = x0{i} + epsilon{i};
        reg = regul_x(epsilon{i}, x0{i},k,i) + reg;
        end
    end
    
    b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
    [~, y_x] = Phi(B_tmp, b, M);
    [~, x_] = Phit(Bt_tmp, y_x, 1);
    [~, dirty_x] = Phit(Bt_tmp, y_x-y, 1);
    [~,phit_y] = Phit(Bt_tmp,y,1);

    
    % Display
     time_tot(nIter+1) = toc;
    data_fid = lambda_scale*sum(abs(y_x(:) - y(:)).^2);
     err_temp = nuo*sum(abs(U1(:) - U2(:)).^2)/2;
    objective(nIter+1) = data_fid + sum(err_temp)/2 + reg; %regul_x(epsilon);
    % keyboard
%     error_dirty(nIter+1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
%  snr_dirty(nIter+1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
%     snr_x(nIter+1) = SNR(x0 + epsilon, x_th);
%     
for i =1:3
    error_dirty(nIter+1,i) = sqrt(sum(sum(dirty_x{i}(:).^2)));
    snr_dirty(nIter+1,i) = SNR(x{i}, phit_y{i});
    snr_x(nIter+1,i) = SNR(x{i}, x_th.x{i});
%      x_reg = x0 + epsilon;
    x_reg{i} = x{i}/max(x{i}(:));
    scale_x(i) = sum(x_reg{i}(:).*x_th.x{i}(:))/sum(x_reg{i}(:).^2);
    scaled_snr(i) = SNR(scale_x(i)*x_reg{i}, x_th.x{i});
end


 figure(100)
     subplot 331, imagesc(x_th.im_true{1}), title('True I'), colorbar, axis image, axis off
     subplot 332, imagesc(x_th.im_true{2}), title('True Q'), colorbar, axis image, axis off  
     subplot 333, imagesc(x_th.im_true{3}), title('True U'), colorbar, axis image, axis off
     subplot 334, imagesc(x{1}), title('Reconstructed I'), colorbar, axis image, axis off
     subplot 335, imagesc(x{2}), title('Reconstructed Q'), colorbar, axis image, axis off
     subplot 336, imagesc(x{3}), title('Reconstructed U'), colorbar, axis image, axis off
     subplot 337, plot(snr_x(:,1))
     subplot 338, plot(snr_x(:,2))
     subplot 339, plot(snr_x(:,3))

snr_x1 = snr_x(nIter+1,2);
%     x_reg = x0 + epsilon;
%     x_reg = x_reg/max(x_reg(:));
%     scale_x = sum(x_reg(:).*x_th(:))/sum(x_reg(:).^2);
%     scaled_snr = SNR(scale_x*x_reg, x_th);
%     
    fprintf('%2s  %2s  %2s\t%2s\t%2s\t%2s  %2s  %2s\n', 'Iteration', 'Objective', 'SNR_I', 'SNR_Q', 'SNR_U', 'scaled SNR', 'Time')
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('%2i  %.2e  %.2e\t%.2e\t%.2e\t%2e  %2e  %2e\n', nIter, objective(nIter+1), snr_x(nIter+1,1), snr_x(nIter+1,2), snr_x(nIter+1,3), scaled_snr, time_tot(nIter+1));
 
fprintf('================================================================================\n');
   
%     if mod(nIter, 5)==0
%         save(['results/img_dde_no_reg_', name], 'U1', 'U2', 'epsilon')
%     end
%     
  %  Global stopping criterion
   if abs((objective(nIter+1) - objective(nIter))/objective(nIter) ) < param_algo.tol_crit
         fprintf('End of algorithm: global stopping criterion satisfied.\n')
         break
     end
%     
   
%     err1 = nuo*sum(abs(U1(:) - U2(:)).^2)/2;
%     err2 = mu*(sum(abs(U1(:) - Phi(:)).^2) + sum(abs(U2(:) - Phi(:)).^2))/2; 
%     data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
%     objective(nIter+1) = data_fid + err1 + err2 + regul_x(epsilon);
  
end

end
