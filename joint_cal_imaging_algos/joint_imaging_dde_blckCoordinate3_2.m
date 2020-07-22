function [U1, U2, epsilon, objective, time_tot, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate3_2(y, Y, x0, epsilon, x_th, u, v, Omega, na, T, J, K, N, S, U1, U2, scale, V, param_im, param_dde, param_algo, W, name, mu,L,Lt,M,param_pd, A_true, At_true,D_th,G_true) %, U1t_,U2t_)


%% Notes:
% J: size of the NFFT support
% S: size of the DDE support
% Ni: image dimension
% Nt: Fourier space dimension (spatial)
% Omega: {T,1}: cell containing the u, v coordinates with redundancy (as for Y matrices, Omega{t}: [na(t), na(t), 2], first slice corresponding to v)
% V: cell containing the gridding coefficients for NFFT, presents redundancy as for Y
% scale: scaling factors associated to the NFFT
% U1 / U2: cell containing the nonzero coefficients of the DDE

%% Initialization

time_tot = zeros(param_dde.max_it+1, 1);
objective = zeros(param_dde.max_it+1, 1);
error_dirty = zeros(param_dde.max_it+1, 1);
snr_dirty = zeros(param_dde.max_it+1, 1);
snr_x = zeros(param_dde.max_it+1, 1);

param_l1.min_x = param_im.min_x;
param_l1.max_x = param_im.max_x;
param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
for i =1:3
    param_l1.mask_app{i} = (abs(x0{i}) > 0);
    param_l1.xA{i} = x0{i};
    if param_pd.dual_fb == 0
        for k=1:param_pd.Ps
            Psitw{k} = param_l1.Psit{k};
            Psiw{k} = param_l1.Psi{k};
            x_th.psit_xo{i,k} = Psitw{k}(x_th.xo{i});
        end
    end
end

param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.weights = 1;
param_l1.eta_o = [1e-2, 1e-1, 1e+5];


lambda_scale = 1/prod(K);
param_algo.eta = param_algo.eta*lambda_scale;
eta = param_algo.eta; 
proxEpsilon = @(x,Ax,i) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1,i) ;
regul_x = @(x,x0,k) param_algo.eta * sum( abs(param_l1.Psit{k}(x0+x)) ) ;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

theta_maxoR = param_dde.theta_maxR;
theta_minoR = param_dde.theta_minR;
theta_maxoI = param_dde.theta_maxI;
theta_minoI = param_dde.theta_minI;
JU2o = param_dde.JU2;
JU1o = param_dde.JU1;
nuo = param_dde.nu;
stop_crit_temp = zeros(T,1);
err_temp = zeros(T,1);
S2 = S^2;
Phi_ = sparse((floor(S2/2)+1)*ones(na(1), 1), (1:na(1)).', ones(na(1), 1), S2, na(1));
Jeps = param_algo.Jeps;
tol_x = param_algo.tol_x;

% Initial value of the objective function
% G = cell(T, 1);

tic

% %par
parfor t = 1:T 
    [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
     err_temp(t) = nuo*sum((abs(U1{t}(:) - U2{t}(:))).^2);
end
toc(tic)

G0_{1} = cell2mat(G1_');
G0_{2} = cell2mat(G2_');
G0_{3} = cell2mat(G3_');
G0_{4} = cell2mat(G4_');

for i = 1:4
% G0{i} = [G{i,1:end}];
B_tmp_{i} = @(x) G0_{i}*so_fft2(x, K, scale);
Bt_tmp{i} = @(x) real(so_fft2_adj((G0_{i})'*x, N, K, scale));
end

for i =1:3
    x{i} = x0{i} + epsilon{i};
end

   b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
%   b = x_th.bright_im;
% [B_x] = Phi(B_tmp_, b, M);

% for i =1:4
% data_fid(i) = lambda_scale*sum(abs(B_x{i} - y(:,i)).^2);
% objective(1) = data_fid + sum(err_temp)/2 + regul_x(epsilon);
% error_dirty(1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
% snr_dirty(1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
% snr_x(i,1) = SNR(b{i}, x_th);
% end

fprintf('================================================================================\n');
fprintf('Joint calibration and imaging (no temporal regularization) [%s]\n', name);
fprintf('================================================================================\n');
%% Global algorithm
tic
for nIter = 1:param_dde.max_it
    
    %%% Update DDEs
    
    for i=1:4
        x_hat{i} = reshape(fftshift(reshape(so_fft2(b{i}, K, scale), K)), [prod(K), 1]); % FT of brightness matrix
    end

     for nIter_dde = 1:param_dde.JUtot
         fprintf('--DDE iterations: %4i / %4i \n', nIter_dde, param_dde.JUtot);
         % Update U1 & U2
         %         par
       parfor t = 1:T % see parfeval or parfor for the update of U1 and U2: compute indices, then send appropriate portion of x_hat
             
             U_old = (U1{t}(:,:,:));
             id = true(na(t),1);
             Ht = cell(8,na(t));
             for a = 1:na(t)
                 id(a) = false;
                 om1 = squeeze(Omega{t}(a,id,1)).';
                 om1(isnan(om1)) = [];
                 om2 = squeeze(Omega{t}(a,id,2)).';
                 om2(isnan(om2)) = [];
                 id1 = indices4Xhat(S, J, K, [om1, om2]);
                 id_nnz = find(~isnan(Y{i,t}(a,id)));
                 [Ht{1,a},Ht{2,a},Ht{3,a},Ht{4,a},Ht{5,a},Ht{6,a},Ht{7,a},Ht{8,a}] = createH1a_32(x_hat, squeeze(V{t}(:, a, id)), U2{t}, J, S, id_nnz,id1,id,size(U2{t}(i,:,id)));
                 id(a) = true;
             end
             
             U1{t} = updateUt122(Y,t, U1{t}, U2{t}, Ht, na(t), JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
             stop_crit_temp(t) = norm(U1{t}(:) - U_old(:))/norm((U_old(:)));
             % %
             %         % update U2{t}
             U_old = U2{t}(:,:,:);
             id = true(na(t),1);
             Ht = cell(8,na(t));
             for a = 1:na(t)
                 id(a) = false;
                 om1 = squeeze(Omega{t}(id,a,1));
                 om1(isnan(om1)) = [];
                 om2 = squeeze(Omega{t}(id,a,2));
                 om2(isnan(om2)) = [];
                 id1 = indices4Xhat(S, J, K, [om1, om2]); % [S2^2*J2*(na-1), na] % /!\ nonzeros systematically returns a column vector
                 id_nnz = find(~isnan(Y{i,t}(id,a)));
                 [Ht{1,a},Ht{2,a},Ht{3,a},Ht{4,a},Ht{5,a},Ht{6,a},Ht{7,a},Ht{8,a}]  = createH2a_32(x_hat, squeeze(V{t}(:, id, a)), U1{t}, J, S, id_nnz,id1,id, size(U1{t}(i,:,id)),id1);
                 id(a) = true;
             end
             U2{t}= updateUt222(Y, t, U1{t}, U2{t}, Ht, na(t), JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
             stop_crit_temp(t) = max(stop_crit_temp(t), norm(U2{t}(:) - U_old(:))/norm(U_old(:)));
             
             % regularization term
             err_temp(t) = nuo*sum((abs(U1{t}(:) - U2{t}(:))).^2);
             %             err_temp(t) = nuo*norm(U1{t}(:) - U2{t}(:)) + mu*(norm(U1{t(:)} - Phi_(:)) + norm(U2{t}(:) - Phi_(:)));
         end
         % DDE stopping criterion
         stop_crit = max(stop_crit_temp) ;
         if stop_crit(end) < param_dde.tol_norm
             fprintf('--DDE iterations: stopping criterion satisfied');
             break
         end
     end
     
    fprintf('--------------------------------------------------------------------------------\n');
    
    %% Update Image
    
    % Create convolution matrix G
        G = cell(T,1);
        % tic
        parfor t = 1:T % parfor
            % id_upper = getIndicesUp(na(t)); % indices corresponding to the upper half of the square matrix V{t}
            [G1_{t},G2_{t},G3_{t},G4_{t}] = createGnufft_T2_parallel_T(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format

        end
        % toc
        G0_{1} = cell2mat(G1_');
        G0_{2} = cell2mat(G2_');
        G0_{3} = cell2mat(G3_');
        G0_{4} = cell2mat(G4_');

        for i = 1:4
            B_tmp{i} = @(x) G0_{i}*so_fft2(x, K, scale);
            Bt_tmp{i} = @(x) real(so_fft2_adj((G0_{i})'*x, N, K, scale));
        end

        for i =1:3
            x{i} = x0{i} + epsilon{i};
        end

%         B_tmp = A_true;
%         Bt_tmp = At_true;

        b_eps = conv_stokes_to_bright([epsilon{1}(:),epsilon{2}(:),epsilon{3}(:)],L,size(epsilon{1}));
        [B_eps,~] = Phi(B_tmp, b_eps, M);
        
        b_xo = conv_stokes_to_bright([x0{1}(:),x0{2}(:),x0{3}(:)],L,size(x0{1}));
        [~,Y_xo] = Phi(B_tmp, b_xo, M);

%%%%%%%%%%

    % Update epsilon
    if nIter == 1
        param_dual.dual_vv = 0;
    end
    
    [epsilon, param_dual] = updateEpsilon2_stokes(y, Y_xo, epsilon, Jeps, tol_x, nIter, proxEpsilon, B_tmp, Bt_tmp, lambda_scale,L, Lt, M, eta, param_l1,param_pd,x_th, param_dual);
    reg = 0;
    
    for i =1:3
        for k =1:param_pd.Ps
        x{i} = x0{i} + epsilon{i};
        reg = regul_x(epsilon{i}, x0{i},k) + reg;
        end
    end
    
    b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
    [~, y_x] = Phi(B_tmp, b, M);
    [~, x_] = Phit(Bt_tmp, y_x, 1);
    [~, dirty_x] = Phit(Bt_tmp, y_x-y, 1);
    [~,phit_y] = Phit(Bt_tmp,y,1);

    % Display monitoring results
    time_tot(nIter+1) = toc;
    data_fid = lambda_scale*sum(abs(y_x(:) - y(:)).^2);
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
 end

end
