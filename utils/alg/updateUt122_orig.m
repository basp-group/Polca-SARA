function Ut1 = updateUt122(Yt, Ut1, Ut2, H1t, na_t, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% U2t: [S2, na]
% U1t: [S2, na]
% Yt: [na, na]
% H1t: [na-1, S2, na]
id = true(na_t,1);
S2 = size(Ut1, 1);
phi = sparse(floor(S2/2)+1, 1, 1, S2, 1);

for a = 1:na_t
    id(a) = false;
    
    % construction of data vector
    Yt_alpha = nonzeros(Yt(a,id));  % transpose of the rows of Y /!\ nonzeros systematically returns a column vector
    u1t = flipud(Ut1(:,a));
    u2t = flipud(Ut2(:,a));
    
    % step size
    Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) H1t(:, :, a)*x, @(x) H1t(:, :, a)'*x, size(u1t)); % lambda_scale*
    gamma1 = 1.9/Lips_temp ;

    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU1o
        % gradient step
        grad1 = 2*lambda_scale*H1t(:, :, a)'*(H1t(:, :, a)*u1t - Yt_alpha);
        grad2 = nuo*(u1t - u2t) + mu*(u1t - phi);
        grad = grad1 + grad2;
        g = u1t - gamma1*grad;
        % proximity step
        vr = min(max(real(g), theta_minoR), theta_maxoR);
        vi = min(max(imag(g), theta_minoI), theta_maxoI);
        u1t = vr + 1i*vi;
    end
    % ----------------------------------------------------
    Ut1(:,a) = flipud(u1t);
    id(a) = true;
end

end

