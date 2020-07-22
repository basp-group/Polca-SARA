function Ut2 = updateUt222_orig(Yt, Ut1, Ut2, Ht2, na_t, JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Ut2: [S2, na]
% Ut1: [S2, na]
% Yt: cell(na,1)
% Ht2: [na-1, S2, na]
id = true(na_t,1);
S2 = size(Ut1, 1);
phi = sparse(floor(S2/2)+1, 1, 1, S2, 1);

for a = 1:na_t
    id(a) = false;
    
    % construction of data vector
    Yt_alpha = nonzeros(Yt(id,a));
    ut1 = conj(Ut1(:,a));
    ut2 = conj(Ut2(:,a));
    
    % step size
    Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) Ht2(:, :, a)*x, @(x) Ht2(:, :, a)'*x, size(ut2)); % lambda_scale*
    gamma2 = 1.9/Lips_temp;
    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU2o
        % gradient step
        grad1 = 2*lambda_scale*Ht2(:, :, a)'*(Ht2(:, :, a)*ut2 - Yt_alpha);
        grad2 = nuo*(ut2 - ut1) + mu*(ut2 - phi);
        grad = grad1 + grad2;
        g = ut2 - gamma2*grad;
        % proximity step
        vr = min(max(real(g), theta_minoR), theta_maxoR);
        vi = min(max(imag(g), theta_minoI), theta_maxoI);
        ut2 = vr + 1i*vi;
    end
    % ----------------------------------------------------
    Ut2(:,a) = conj(ut2);
    id(a) = true;
end

end
