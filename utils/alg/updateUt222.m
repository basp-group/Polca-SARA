function Ut2__ = updateUt222(Yt_, t, Ut1_, Ut2_, H2t, na_t, JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here

% Ut2: [S2, na]
% Ut1: [S2, na]
% Yt: cell(na,1)
% Ht2: [na-1, S2, na]
id = true(na_t,1);
S2 = size(Ut1_, 1);
phi = sparse(floor(S2/2)+1, 1, 1, S2, 1);
z1 = size(Ut1_(1,:,:));
z2 = size(Ut2_(1,:,:));

% Function handles to create H matrix for each measurement made by the
% correlator
H{1} = @(u1,u2,a) (H2t{1,a}+ H2t{2,a})*u1 + (H2t{3,a}+ H2t{4,a})*u2;
H{3} = @(u1,u2,a) (H2t{5,a}+ H2t{6,a})*u1 + (H2t{7,a}+ H2t{8,a})*u2;
H{2} = H{1};
H{4} = H{3};

% Adjoint operators
Ht{1} = @(y,a) [(H2t{1,a}+H2t{2,a})'*y, (H2t{3,a}+H2t{4,a})'*y].';
Ht{3} = @(y,a) [(H2t{5,a}+H2t{6,a})'*y, (H2t{7,a}+H2t{8,a})'*y].';
Ht{2} = Ht{1};
Ht{4} = Ht{3};

for i = 1:4
    Yt{i} = Yt_{i,t};
    Ut1{i} = reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ut2{i} = reshape(Ut2_(i,:,:),[z2(2:end) 1]);
end

for a = 1:na_t
    id(a) = false;
    for i = 1:4
        
        % construction of data vector
        Yt_alpha{i} = nonzeros(Yt{i}(id,a));
        u1t{i} = conj(Ut1{i}(:,a));
        u2t{i} = conj(Ut2{i}(:,a));
        % step size
        Lips_temp(i) = mu + nuo + 2*lambda_scale*pow_method_stokes_cal(@(x1,x2) H{i}(x1,x2,a), @(x) Ht{i}(x,a), size(u2t{i})); % lambda_scale*
        gamma2(i) = 1.9/Lips_temp(i) ;
        
    end
    
    %     Lips_temp = mu + nuo + 2*lambda_scale*pow_method(@(x) Ht2(:, :, a)*x, @(x) Ht2(:, :, a)'*x, size(ut2)); % lambda_scale*
    %     gamma2 = 1.9/Lips_temp;
    
    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU2o
        % gradient step
        %         grad1 = 2*lambda_scale*Ht2(:, :, a)'*(Ht2(:, :, a)*ut2 - Yt_alpha);
        grad1{1} = 2*lambda_scale*Ht{1}((H{1}(u2t{1},u2t{3},a) - Yt_alpha{1}),a); %u2t{2}
        grad1{2} = 2*lambda_scale*Ht{2}((H{2}(u2t{2},u2t{4},a) - Yt_alpha{2}),a); %u2t{3}
        grad1{3} = 2*lambda_scale*Ht{3}((H{3}(u2t{1},u2t{3},a) - Yt_alpha{3}),a); %u2t{2}
        grad1{4} = 2*lambda_scale*Ht{4}((H{4}(u2t{3},u2t{4},a) - Yt_alpha{4}),a); 
        
        gu{1} = grad1{1}(1,:) + grad1{3}(1,:);
        gu{2} = grad1{1}(2,:) + grad1{3}(2,:);
        gu{3} = grad1{2}(1,:) + grad1{4}(1,:);
        gu{4} = grad1{2}(2,:) + grad1{4}(2,:);
        
        
        for i = 1:4
            grad2{i} = nuo*(u2t{i} - u1t{i}); % + mu*(u2t{i} - phi{i});
            grad{i} = gu{i}.' + grad2{i};
            g{i} = u2t{i} - gamma2(i)*grad{i};
            % proximity step
            vr{i} = min(max(real(g{i}), theta_minoR(:,i)), theta_maxoR(:,i));
            vi{i} = min(max(imag(g{i}), theta_minoI(:,i)), theta_maxoI(:,i));
            u2t{i} = vr{i} + 1i*vi{i};
        end
        
    end
    % ----------------------------------------------------
    for i =1:4
%         Ut2__ = u2t{i}; 
         Ut2__(i,:,a) = conj(u2t{i});
    end
    id(a) = true;
end

end
