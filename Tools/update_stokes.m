function [x, snr] = update_stokes(y, Jeps, tol_x, n_it,  B_tmp, Bt_tmp, L, Lt, M, eta, param_l1,param_pd,x_th, param_dual, Ax, epsilon,eta_)
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

SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

param_pd.evl = Ax;
param_pd.nu3 = Ax;
data_fid = 0;
l1_term = 0;
reweight_step_count = 1;
step = 1500;
q_rw = 500;

[~,grad_x] = Phit(Bt_tmp,y,2,L,Lt);
for i = 1:3
    x{i} = epsilon{i}; %grad_x{i}/Ax;
    x{i} = reshape(x{i},size(x_th.x{i}));
%       x{i} = x_th.x{i};
end

fprintf('================================================================\n');
fprintf('Imaging with approximate DDEs \n');
fprintf('================================================================\n');

profile on

% tm = tic;
for q = 1:Jeps % around 4 s/iteration
    x_old = x;
     
    b_x= conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1}));
    [~,y_x] = Phi(B_tmp, b_x, M);
   
   [~,grad_x] = Phit(Bt_tmp,y_x-y,2,L,Lt);
   
   data_fid(q) = 0.5*norm(y(:)-y_x(:));

    if param_pd.dual_fb
            l1_term(q) = eta(1)*sum(param_l1.weights{1}(:).*abs(param_l1.Psit(x{1}))) + eta(2)*sum(param_l1.weights{2}(:).*abs(param_l1.Psit(x{2}))) + eta(3)*sum(param_l1.weights{3}(:).*abs(param_l1.Psit(x{3})));
% tm = tic;

        
%par
for i = 1:3
            x_tmp{i} = x{i} - (1.9/Ax)*real(reshape(real(grad_x{i}),size(x{i}))); % Gradient step
             x{i} = solver_prox_L1_full_image(x_tmp{i}, eta(i)*1.9/Ax, param_l1,i) ; % eta(i)*1.9/Ax
% x{i} = x_tmp{i};
            snr{i}(q) = SNR(x{i}, x_th.x{i});
        end
% %         tm = toc(tm);
% %  fprintf('Time : %3.3f\n\n\n',tm);
    else
        for i = 1:3
           x_tmp{i} =x{i} - (1.9/Ax)*real(reshape(real(grad_x{i}),size(x{i}))); % Gradient step
        end
      
        if q == 1 && n_it == 1
            P = 3; R = 1;
            param = param_pd;
            
            def_var_real_data; % to define and initialize variables
        end
% tm = tic;

        [x, param_dual,l1]  = pdfb_stokes_imaging_cal_ddes(x_tmp, 3, param_l1, param_pd, x_th,param_dual,eta_);
%         tm = toc(tm);
%         fprintf('Time : %3.3f\n\n\n',tm);
        for i = 1:3
            snr{i}(q) = SNR(x{i}, x_th.x{i});
        end
            l1_term(q) = l1;

    end
     crit(q) = data_fid(q) + l1_term(q);

%      toc(tic);
if rem(q,200) == 0
     figure(200)
     subplot 331, imagesc(log10(abs(x_th.x{1}))), title('True I'), colorbar, axis image, axis off, caxis([-3.8 -0.2])
     subplot 332, imagesc(log10(abs(x_th.x{2}))), title('True Q'), colorbar, axis image, axis off,  caxis([-3.8 -0.2])
     subplot 333, imagesc(log10(abs(x_th.x{3}))), title('True U'), colorbar, axis image, axis off, caxis([-3.8 -0.2])
     subplot 334, imagesc(log10(abs(x{1}))), title('Reconstructed I'), colorbar, axis image, axis off, caxis([-3.8 -0.2])
     subplot 335, imagesc(log10(abs(x{2}))), title('Reconstructed Q'), colorbar, axis image, axis off, caxis([-3.8 -0.2])
     subplot 336, imagesc(log10(abs(x{3}))), title('Reconstructed U'), colorbar, axis image, axis off, caxis([-3.8 -0.2])
     subplot 337, plot(snr{1})
     subplot 338, plot(snr{2})
     subplot 339, plot(snr{3})
%      
     figure(400)
     subplot 131, semilogy(data_fid), title('data fidelity')
     subplot 132, semilogy(l1_term), title('l1 term')
     subplot 133, semilogy(crit), title('objective function')
     pause(0.1)
     
end
    % stopping criterion (imaging step)
if (q>q_rw+step) && (norm([x{1:end}] - [x_old{1:end}], 'fro')/norm([x{1:end}], 'fro') < tol_x)
        disp(['x: stopping criterion reached, glob. iteration = ', num2str(n_it)])
        disp(['x: stopping criterion reached, inner iteration = ', num2str(q)])
%         break
%     end
% end
% tm = toc(tm);
%  fprintf('Time : %3.3f\n\n\n',tm);
%  


 % Computation of weights for reweighting scheme
  fprintf('\n\n\n\n\n\n\n Performing reweight no %d \n\n\n\n\n', reweight_step_count);
   reweight_step_count = reweight_step_count+1;
   q_rw = q;
    if param_pd.dual_fb
        for i = 1:3
            d_val = abs(param_l1.Psit(reshape(x{i},size(x_th.x{i}))));
            d_val_max = max(d_val(:));
            if d_val_max == 0
                d_val_max =1;
            end
            param_l1.weights{i} = param_l1.reweight_alpha ./ (param_l1.reweight_alpha + (d_val./d_val_max));
            param_l1.weights{i}(d_val > d_val_max * param_l1.reweight_abs_of_max) = 0;
%             eta{i} = param_l1.eta_o(i)*weights{i}*1.9/Ax;
            param_l1.reweight_alpha = param_l1.reweight_alpha_ff .* param_l1.reweight_alpha;
        end
    else
        for i = 1:3
            for k = 1:param_pd.Ps
                d_val = abs(param_l1.Psit{k}(reshape(x{i},size(x_th.x{i}))));
                d_val_max = max(d_val(:));
                if d_val_max == 0
                    d_val_max =1;
                end
                 param_l1.weights{i,k} = param_l1.reweight_alpha ./ (param_l1.reweight_alpha + (d_val./d_val_max));
                 param_l1.weights{i,k}(d_val > d_val_max * param_l1.reweight_abs_of_max) = 0;
                eta_{i,k} = param_l1.eta_o(i)*param_l1.weights{i,k};
                param_l1.reweight_alpha = param_l1.reweight_alpha_ff .* param_l1.reweight_alpha;
            end
        end
 end
 
end
end
end
