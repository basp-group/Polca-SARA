function [x, snr, crit] = imaging_fb_th(y, x_true, G0, sp_scale, Ax, N, K, nIter, param_algo, A, At, stokes_P, L, Lt,Ni)
%%

% Output metrics
SNR = @(x, xtrue) 10*log10(sum(xtrue(:).^2)/sum((x(:) - xtrue(:)).^2));

%% Image initialization from pre-calibrated DIEs

if ~isfield(param_algo, 'stop_crit'), param_algo.stop_crit = 1e-5; end
% A = @(x) G0*so_fft2(x, K, sp_scale);
% At = @(x) real(so_fft2_adj(G0'*x, N, K, sp_scale));

%%Added by jb
% x = 2*At(y)/Ax; % to be verified
for i = 1:4
    b{i} = 2*At{i,1}(y(:,i))/Ax;
    b{i} = b{i}(:);
    
end

  x = mat2cell([b{1:4}]*Lt,size(b{1},1), [1 1 1])'; % conversion from brightness matrix to Stokes matrix
for i=1:stokes_P
    x{i} = reshape(x{i},Ni(2),Ni(1));
end
  
  b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,Ni);
  
eta_l1 = param_algo.eta;
fprintf('================================================================================\n');
fprintf('Imaging with true DDEs (eta = %5e)\n', eta_l1);
fprintf('================================================================================\n');

% if strcmp(param_algo.opt_dict, 'Dirac')
%     prox_func =@(x) max(proxl1(x, eta_l1*(1.9/Ax) ), 0) ;
% else
%     param_l1.Psi = param_algo.Psi ;
%     param_l1.Psit = param_algo.Psit ;
%     param_l1.real = 1 ;
%     param_l1.pos = 1 ;
%     param_l1.verbose = 0 ;
%     param_l1.weights = 1 ;
%     prox_func =@(x) solver_prox_L1(x, eta_l1*(1.9/Ax), param_l1) ;
% end

param_l1.Psi = param_algo.Psi;
param_l1.Psit = param_algo.Psit;
param_l1.real = 1;
param_l1.pos = 1;
param_l1.verbose = 0;
param_l1.weights = 1;
param_l1.im_size = Ni;
prox_func = @(x,pos) solver_prox_L1(x, eta_l1*(1.9/Ax), param_l1,pos);
% objective = zeros(nIter, 1);
% objective(1) = sum(abs(A(x) - y).^2) + eta_l1*sum(abs(param_algo.Psit(x)));



tic
for q = 1:nIter
    xold = x;
    bold = b;
    
    for i=1:4
    grad_b{i} = 2*At{i,1}(A{i,1}(b{i}) - y(:,i));
    grad_b{i} = grad_b{i}(:);
    end
    
    % Conversion from brightness to Stokes domain
    grad = mat2cell([grad_b{1:4}]*Lt,size(grad_b{1},1), [1 1 1])';
    
    l1_val(q) = 0;
    
    for i=1:stokes_P
        if i>1
            param_l1.pos = 0;
        else
            param_l1.pos = 1;
        end
    x{i} = prox_func(x{i} - (1.9/Ax)*reshape(grad{i},Ni(2),Ni(1)),param_l1.pos); 
    l1_val(q) = eta_l1*sum(abs(param_algo.Psit(x{i})))+ l1_val(q);
    
    crit_stop{i}(q) = norm(x{i}(:) - xold{i}(:))/norm(xold{i}(:));
    if crit_stop{i}(q) < param_algo.stop_crit
        disp(['Stopping criterion reached: q = ', num2str(q)])
        break
    end
    
    snr{i}(q) = SNR(x{i}, x_true{i});
    end
    
    b = conv_stokes_to_bright([x{1}(:),x{2}(:),x{3}(:)],L,size(x{1})); %Conversion from Stokes to brightness matrix
    
    obj_data_val(q) = 0;
    for i =1:4
        obj_data{i} = A{i,1}(b{i}) - y(:,i);
        obj_data_val(q) = sum(abs(obj_data{i}).^2) + obj_data_val(q);
    end
    
%         objective(q+1) = sum(abs(A(x) - y).^2) + eta_l1*sum(abs(param_algo.Psit(x)));    
    objective(q) = obj_data_val(q) + l1_val(q);    
    
    % Figures
    figure(200)
    subplot 331, semilogy(objective)
    subplot 334, semilogy(crit_stop{1})
    subplot 335, semilogy(crit_stop{2})
    subplot 336, semilogy(crit_stop{3})
    subplot 337, plot(snr{1})
    subplot 338, plot(snr{2})
    subplot 339, plot(snr{3})

pause(0.1)
end
time_img = toc;


crit = sum(abs(A(x) - y).^2) + eta_l1*sum(abs(param_algo.Psit(x)));
% error_dirty = sqrt(sum(sum((At(A(x) - y)).^2)));


x1 = x/max(x(:));
scale_x = sum(x1(:).*x_true(:))/sum(x1(:).^2);
scaled_snr = SNR(scale_x*x1, x_true);

% Display results
fprintf('%10s\t%15s\t%15s\t%15s\t%15s\n', 'Iteration', 'Objective', 'SNR', 'scaled SNR', 'Time')
fprintf('--------------------------------------------------------------------------------\n');
fprintf('%10i\t%15e\t%15e\t%15e\t%15e\n', q, crit, snr, scaled_snr, time_img);
fprintf('================================================================================\n');

end
