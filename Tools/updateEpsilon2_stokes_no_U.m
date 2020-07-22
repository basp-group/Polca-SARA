 function [epsilon, param_dual] = updateEpsilon2_stokes_no_U(y, y0, epsilon, Jeps, tol_x, n_it, proxEpsilon, B_tmp, Bt_tmp, lambda_scale, L, Lt, M, eta, param_l1,param_pd,x_th, param_dual)
% Update the image parameter epsilon (PALM descent steps).
%-------------------------------------------------------------------------%
% Input:
% > y            : data vector [M, 1]
% > x0           : reference image term [N]
% > epsilon      : initial epsilon term [N]
% > Jeps         : maximum number of iterations
% > tol_x        : stopping criterion
% > proxEpsilon  : lambda function for the proximal operator
% > B_tmp        : lambda function representing the direct measurement operator
% > Bt_tmp       : adjoint of B_tmp (lambda function)
% > lambda_scale : rescaling term (used everywhere in the problem)
%                  lambda_scale = 1/prod(K), K dimension of the spatial
%                  Fourier space
%
% Output:
% < epsilon   : updated epsilon term
%
%-------------------------------------------------------------------------%

     Ax = 2*lambda_scale*op_norm_stokes_RIME_no_U(B_tmp,Bt_tmp, size(epsilon{1}), 1e-6, 200, 10, 3,1,L,Lt,M); % operator norm, the factor 2 comes from the complex case

%Ax = 2*8.018278e+06;
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

param_pd.evl = Ax;
% param_pd.nu3 = Ax;

ytmp = y - y0;

% epsilon = x_th.im_true;
stokes_P = 2;

tm = tic;
for q = 1:Jeps % around 4 s/iteration
    eps_old = epsilon;
    
    b_eps = conv_stokes_to_bright([epsilon{1}(:),epsilon{2}(:)],L,size(epsilon{1}));
    [~,y_eps] = Phi(B_tmp, b_eps, M);
   
%     [~,grad_eps] = Phit(Bt_tmp,y_eps - ytmp,2*lambda_scale);
      [~,grad_eps] = Phit_no_U(Bt_tmp,y_eps-ytmp,1,L,Lt);
%     grad_eps = 2*lambda_scale*Bt_tmp(B_tmp(epsilon) - ytmp);

    data_fid(q) = 0.5*norm(ytmp(:)-y_eps(:));
%     l1_term(q) = (1.9/Ax)*(eta(1)*sum(abs(param_l1.Psit(epsilon{1}+x_th.x_approx{1}))) + eta(2)*sum(abs(param_l1.Psit(epsilon{2}+x_th.x_approx{2}))) + eta(3)*sum(abs(param_l1.Psit(epsilon{3}+x_th.x_approx{3}))));


    if param_pd.dual_fb == 1
        for i = 1:stokes_P
            eps_tmp{i} = epsilon{i} - (1.9/Ax)*2*lambda_scale*real(reshape(grad_eps{i},size(epsilon{i}))); % Gradient step
            epsilon{i} = solver_prox_L1_full_image_no_U(eps_tmp{i}, eta(i)*1.9/Ax, param_l1,i) ;
            snr{i}(q) = SNR(epsilon{i}, x_th.eps_true{i});
%             snr_{i}(q) = SNR(epsilon{i}+x_th.x_approx{i}, x_th.true_best{i});
%             snr__{i}(q) = SNR(epsilon{i}+x_th.x_approx{i}, x_th.im_true{i});

        end
    else
        for i = 1:stokes_P
            eps_tmp{i} = epsilon{i} - (1.9/Ax)*2*lambda_scale*real(reshape(grad_eps{i},size(epsilon{i}))); % Gradient step
        end
      
        if q == 1 && n_it == 1
            P = 3; R = 1;
            param = param_pd;
            
            def_var_real_data; % to define and initialize variables
        end
%         n_itt = (n_it == 1 && q == 1);

        [epsilon, param_dual]  = pdfb_stokes_imaging_cal_ddes_no_U(eps_tmp, stokes_P, param_l1, param_pd, x_th,param_dual, eta);
        %  (eps_tmp,3, param_l1,param_pd,x_th); % Proximal step using primal-dual
        for i = 1:stokes_P
             snr{i}(q) = SNR(epsilon{i}, x_th.eps_true{i});
%              snr_{i}(q) = SNR(epsilon{i}+x_th.x_approx{i}, x_th.true_best{i});
%             snr__{i}(q) = SNR(epsilon{i}+x_th.x_approx{i}, x_th.im_true{i});
% 
         end
    end
% %      crit(q) = data_fid(q) + l1_term(q);

%      toc(tic);
if ~rem(q,10)
     figure(600)
     subplot 331, imagesc(x_th.eps_true{1}), title('True I'), colorbar, axis image, axis off
     subplot 332, imagesc(x_th.eps_true{2}), title('True Q'), colorbar, axis image, axis off  
     subplot 333, imagesc(x_th.eps_true{3}), title('True U'), colorbar, axis image, axis off
       subplot 334, imagesc(epsilon{1}), title('Reconstructed I'), colorbar, axis image, axis off
     subplot 335, imagesc(epsilon{2}), title('Reconstructed Q'), colorbar, axis image, axis off
%      subplot 336, imagesc(epsilon{3}), title('Reconstructed U'), colorbar, axis image, axis off
     subplot 337, plot(snr{1})
     subplot 338, plot(snr{2})
%       subplot 339, plot(snr{3})
     figure(250)
      subplot 331, semilogy(data_fid), title('data fidelity')
%      subplot 332, semilogy(l1_term), title('l1 term')
%      subplot 333, semilogy(crit), title('objective function')
%      subplot 334, plot(snr_{1})
%      subplot 335, plot(snr_{2})
%      subplot 336, plot(snr_{3})
%      subplot 337, plot(snr__{1})
%      subplot 338, plot(snr__{2})
%      subplot 339, plot(snr__{3})
     pause(0.1)
end
    % stopping criterion (imaging step)
    if (q>10) && (norm([epsilon{1:end}] - [eps_old{1:end}], 'fro')/norm([epsilon{1:end}], 'fro') < tol_x)
        disp(['x: stopping criterion reached, glob. iteration = ', num2str(n_it)])
        disp(['x: stopping criterion reached, inner iteration = ', num2str(q)])
        break
    end
end
tm = toc(tm);
 fprintf('Time : %3.3f\n\n\n',tm);
end
