function val = op_norm_stokes_RIME(A,At,im_size, tol, max_iter,T,stokes_P)
% op_norm_stokes(A, At, T, W, im_size, tol, max_iter, verbose,P,R)
%% computes the maximum eigen value of the compund operator At*A


im_size = [im_size, stokes_P];
x = randn(im_size);
% x_mat = cell2mat(x);
x = x/norm(x(:));

% x = x/norm(x(:));
init_val = 1;

for k = 1:max_iter
%     for i = 1:P
%         xx{i} = zeros(size(W{1}{1},1),1);
%         ns{i} = A(x(:,:,i));
%         for q=1:R
%             y{i}{q} = T{i}{q}*ns{i}(W{i}{q});
%             
%             Tt{i,q} = T{i}{q}';
%           
%             xx{i}(W{i}{q}) = xx{i}(W{i}{q}) + Tt{i,q}*y{i}{q};
%         end
%         
%             x(:,:,i) = real(At(xx{i}));
%         end
%  

for i = 1:stokes_P
    for t = 1:T
        
    x(:,:,i)=At{i,t}(A{i,t}(x(:,:,i)));
    
    end
end

val = norm(x(:));
%     val = norm(x(:));
    rel_var = abs(val-init_val)/init_val;
    if (verbose > 1)
        fprintf('Iter = %i, norm = %e \n',k,val);
    end
    if (rel_var < tol)
       break;
    end
    init_val = val;
    x = x/val;
    
end

if (verbose > 0)
    fprintf('Norm = %e \n\n', val);
end

end
