function epsilon = updateEpsilon2(y, x0, epsilon, Jeps, tol_x, nit_tot, proxEpsilon, B_tmp, Bt_tmp, lambda_scale)
%UNTITLED10 Summary of this function goes here
%   Detailed explanation goes here

Ax = 2*lambda_scale*pow_method(B_tmp, Bt_tmp, size(epsilon));
ytmp = y - B_tmp(x0);

for q = 1:Jeps % around 4 s/iteration
    eps_old = epsilon;
    grad_eps = 2*lambda_scale*Bt_tmp(B_tmp(epsilon) - ytmp);
    eps_tmp = epsilon - (1.9/Ax)*real(grad_eps);
    epsilon = proxEpsilon(eps_tmp, Ax);
     
    % stopping criterion (imaging step)
    if (q>10) && (norm(epsilon - eps_old, 'fro')/norm(epsilon, 'fro') < tol_x)
        disp(['x: stopping criterion reached, glob. iteration = ', num2str(nit_tot)])
        disp(['x: stopping criterion reached, inner iteration = ', num2str(q)])
        break
    end
end

end
