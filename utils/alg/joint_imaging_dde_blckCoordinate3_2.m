function [U1, U2, epsilon, objective, time_tot, error_dirty, snr_dirty, snr_x] = joint_imaging_dde_blckCoordinate3_2(y, Y, x0, epsilon, x_th, u, v, Omega, na, T, J, K, N, S, U1, U2, scale, V, param_im, param_dde, param_algo, W, name, mu)


%% Remarks:
% - remove as many structures as possible (avoid unnecessary communication when distributing the computation)
% - data generation to be modified (to be seen)
% - parfeval to be investigated
% - make the interface with Arwa's main file (data generation and formatting, ...)
% - add code to implement the temporal prior (temporal NUFFT)

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
param_l1.mask_app = (x0 > 0);
param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.weights = 1;
param_l1.xA = x0;

lambda_scale = 1/prod(K);
param_algo.eta = param_algo.eta*lambda_scale;
proxEpsilon = @(x,Ax) solver_prox_L1_full_image(x, param_algo.eta*1.9/Ax, param_l1) ;
regul_x = @(x) param_algo.eta * sum( abs(param_l1.Psit(x0+x)) ) ;
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
Phi = sparse((floor(S2/2)+1)*ones(na(1), 1), (1:na(1)).', ones(na(1), 1), S2, na(1));

% Initial value of the objective function
G = cell(T, 1);
parfor t = 1:T % parfor
%     id_upper = getIndicesUp(na(t)); % indices corresponding to the upper half of the square matrix V{t}
    G{t} = createGnufft_T2(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:, id_upper) see st.uu..., see format
    % err_temp(t) = nuo*sum((abs(U1{t}(:) - U2{t}(:))).^2);
    err_temp(t) = nuo*sum(abs(U1{t}(:) - U2{t}(:)).^2) + mu*(sum(abs(U1{t}(:) - Phi(:)).^2) + sum(abs(U2{t}(:) - Phi(:)).^2));
end
G0 = cell2mat(G);
B_tmp = @(x) G0*so_fft2(x, K, scale);
Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));


data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
objective(1) = data_fid + sum(err_temp)/2 + regul_x(epsilon);
error_dirty(1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
snr_dirty(1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
snr_x(1) = SNR(x0 + epsilon, x_th);

fprintf('================================================================================\n');
fprintf('Joint calibration and imaging (no temporal regularization) [%s]\n', name);
fprintf('================================================================================\n');
%% Global algorithm
tic
for nIter = 1:param_dde.max_it
    
    %%% Update DDEs
    % x_hat = reshape(fftshift(fft2((x0 + epsilon).*scale, K(1), K(2))), [prod(K), 1]); % 2D format, might be modified later, 0 in the center  
    x_hat = reshape(fftshift(reshape(so_fft2(x0 + epsilon, K, scale), K)), [prod(K), 1]);
    
    for nIter_dde = 1:param_dde.JUtot
        fprintf('--DDE iterations: %4i / %4i \n', nIter_dde, param_dde.JUtot);
        % Update U1 & U2
        parfor t = 1:T % see parfeval or parfor for the update of U1 and U2: compute indices, then send appropriate portion of x_hat
            % tic
            U_old = U1{t};          
            id = true(na(t),1);
            Ht = zeros(na(t)-1, S2, na(t));
            for a = 1:na(t)
                id(a) = false;
                % id1 = indices4Xhat(S, J, Nt, [nonzeros(squeeze(Omega{t}(a,id,1))), nonzeros(squeeze(Omega{t}(a,id,2)))]); % [S2^2*J2*(na-1), na] % /!\ nonzeros systematically returns a column vector
                % id_nnz = find(Y{t}(a,id));
                om1 = squeeze(Omega{t}(a,id,1)).';
                om1(isnan(om1)) = [];
                om2 = squeeze(Omega{t}(a,id,2)).';
                om2(isnan(om2)) = [];
                id1 = indices4Xhat(S, J, K, [om1, om2]);
                id_nnz = find(~isnan(Y{t}(a,id)));
                Ht(:, :, a) = createH1a_32(x_hat(id1), squeeze(V{t}(:, a, id)), U2{t}(:,id), J, S, id_nnz).';         
                id(a) = true;
            end 
            % U1{t} = updateUt121(Y{t}, U1{t}, U2{t}, Ht, na(t), JU1o, nuo, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            U1{t} = updateUt122(Y{t}, U1{t}, U2{t}, Ht, na(t), JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            stop_crit_temp(t) = norm(U1{t} - U_old, 'fro')/norm(U_old, 'fro');
            
            % update U2{t}
            U_old = U2{t};
            id = true(na(t),1);
            Ht = zeros(na(t)-1, S2, na(t));
            for a = 1:na(t)
                id(a) = false;
                %id1 = indices4Xhat(S, J, Nt, [nonzeros(squeeze(Omega{t}(id,a,1))), nonzeros(squeeze(Omega{t}(id,a,2)))]); % [S2^2*J2*(na-1), na] % /!\ nonzeros systematically returns a column vector
                %id_nnz = find(Y{t}(id,a));
                om1 = squeeze(Omega{t}(id,a,1));
                om1(isnan(om1)) = [];
                om2 = squeeze(Omega{t}(id,a,2));
                om2(isnan(om2)) = [];
                id1 = indices4Xhat(S, J, K, [om1, om2]); % [S2^2*J2*(na-1), na] % /!\ nonzeros systematically returns a column vector
                id_nnz = find(~isnan(Y{t}(id,a)));
                Ht(:, :, a) = createH2a_32(x_hat(id1), squeeze(V{t}(:, id, a)), U1{t}(:,id), J, S, id_nnz);  
                id(a) = true;
            end 
            %U2{t} = updateUt221(Y{t}, U1{t}, U2{t}, Ht, na(t), JU2o, nuo, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            U2{t} = updateUt222(Y{t}, U1{t}, U2{t}, Ht, na(t), JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI);
            stop_crit_temp(t) = max(stop_crit_temp(t), norm(U2{t} - U_old, 'fro')/norm(U_old, 'fro'));
            
            % regularization term
            % err_temp(t) = nuo*sum((abs(U1{t}(:) - U2{t}(:))).^2);
            err_temp(t) = nuo*sum(abs(U1{t}(:) - U2{t}(:)).^2) + mu*(sum(abs(U1{t}(:) - Phi(:)).^2) + sum(abs(U2{t}(:) - Phi(:)).^2));
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
        G{t} = createGnufft_T2(U1{t}, U2{t}, [v{t}, u{t}], K, S, J, W{t}); % V{t}(:,id_upper)
    end
    % toc
    G0 = sparse(cell2mat(G));

    % Update epsilon
    B_tmp = @(x) G0*so_fft2(x, K, scale);
    Bt_tmp = @(x) real(so_fft2_adj(G0'*x, N, K, scale));
    epsilon = updateEpsilon2(y, x0, epsilon, param_algo.Jeps, param_algo.tol_x, nIter, proxEpsilon, B_tmp, Bt_tmp, lambda_scale,x_th);
    
    % Display monitoring results
    time_tot(nIter+1) = toc;
    data_fid = lambda_scale*sum(abs(B_tmp(x0 + epsilon) - y).^2);
    objective(nIter+1) = data_fid + sum(err_temp)/2 + regul_x(epsilon);
    % keyboard
    error_dirty(nIter+1) = sqrt(sum(sum((Bt_tmp(B_tmp(x0 + epsilon) - y)).^2)));
    snr_dirty(nIter+1) = SNR(Bt_tmp(B_tmp(x0 + epsilon)), Bt_tmp(y));
    snr_x(nIter+1) = SNR(x0 + epsilon, x_th);
    
    x_reg = x0 + epsilon;
    x_reg = x_reg/max(x_reg(:));
    scale_x = sum(x_reg(:).*x_th(:))/sum(x_reg(:).^2);
    scaled_snr = SNR(scale_x*x_reg, x_th);
    
    fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n', 'Iteration', 'Objective', 'SNR', 'scaled SNR', 'Time')
    fprintf('--------------------------------------------------------------------------------\n');
    fprintf('%10i\t%15e\t%15e\t%15e\t%15e\n', nIter, objective(nIter+1), snr_x(nIter+1), scaled_snr, time_tot(nIter+1));
    fprintf('================================================================================\n');
   
    if mod(nIter, 5)==0
        save(['results/img_dde_no_reg_', name], 'U1', 'U2', 'epsilon')
    end
    
    % Global stopping criterion
    if abs((objective(nIter+1) - objective(nIter))/objective(nIter) ) < param_algo.tol_crit
        fprintf('End of algorithm: global stopping criterion satisfied.\n')
        break
    end
end

end
