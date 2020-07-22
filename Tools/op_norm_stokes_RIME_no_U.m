function val = op_norm_stokes_RIME_no_U(A,At,im_size, tol, max_iter,T,stokes_P,verbose,L,Lt,M)
% op_norm_stokes(A, At, T, W, im_size, tol, max_iter, verbose,P,R)
%% computes the maximum eigen value of the compund operator At*A


im_size = [im_size, stokes_P];
x = randn(im_size);
% x_mat = cell2mat(x);
x = x/norm(x(:));
% x = x/norm(x(:));
init_val = 1;

for k = 1:max_iter
    if k==1
 x1 = x(:,:,1);
 x2 = x(:,:,2);
%  x3 = x(:,:,3);
    else
        x1 = x(:,1);
        x2 = x(:,2);
%         x3 = x(:,3);
    end
    St = [x1(:),x2(:)];
    B1 = St*L; % Brightness matrix

 for i =1:4
    B{i} = reshape(B1(:,i), im_size(2), im_size(1)); 
 end
    
    y_ = Phi(A,B,M);
    
    for i = 1:4
        yy(:,i) = y_{i};
    end
    x_ = Phit_no_U(At,yy,1,L,Lt);
%     x_{i}= At{i}(A{i}(B{i}));
%     x_{i} = x_{i}(:);
   
%  end
    x = [x_{1:4}]*Lt;
    x = real(x);
%     x = ([x_{1:4}]*Lt,size(x_{1},1), [1 1 1])';


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



