function Ut1__ = updateUt122(Yt_, t, Ut1_, Ut2_, H1t, na_t, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% U2t: [S2, na]
% U1t: [S2, na]
% Yt: [na, na]
% H1t: [na-1, S2, na]
id = true(na_t,1);
S2 = size(Ut1_, 2);
phi = sparse(floor(S2/2)+1, 1, 1, S2, 1);
z1 = size(Ut1_(1,:,:));
z2 = size(Ut2_(1,:,:));

% Function handles to create H matrix for each measurement made by the
% correlator
H{1} = @(u1,u2,a) (H1t{1,a}+ H1t{3,a})*u1 + (H1t{2,a}+ H1t{4,a})*u2;
H{2} = @(u1,u2,a) (H1t{5,a}+ H1t{7,a})*u1 + (H1t{6,a}+ H1t{8,a})*u2;
H{3} = H{1};
H{4} = H{2};

% Adjoint operators
Ht{1} = @(y,a) [(H1t{1,a}+H1t{3,a})'*y, (H1t{2,a}+H1t{4,a})'*y].';
Ht{2} = @(y,a) [(H1t{5,a}+H1t{7,a})'*y, (H1t{6,a}+H1t{8,a})'*y].';
Ht{3} = Ht{1};
Ht{4} = Ht{2};

for i = 1:4
    Yt{i} = Yt_{i,t};
    Ut1{i} = reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ut2{i} = reshape(Ut2_(i,:,:),[z2(2:end) 1]);
end

for a = 1:na_t
    id(a) = false;
    
    for i = 1:4
    % construction of data vector
    Yt_alpha{i} = nonzeros(Yt{i}(a,id));  % transpose of the rows of Y /!\ nonzeros systematically returns a column vector
    u1t{i} = flipud(Ut1{i}(:,a));
    u2t{i} = flipud(Ut2{i}(:,a));
    Lips_temp(i) = nuo + 2*lambda_scale*pow_method_stokes_cal(@(x1,x2) H{i}(x1,x2,a), @(x) Ht{i}(x,a), size(u1t{i})); % + muo 
    gamma1(i) = 1.9/Lips_temp(i) ;
    end
  
    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU1o
        % gradient step
%         grad1 = 2*lambda_scale*H1t(:, :, a)'*(H1t(:, :, a)*u1t - Yt_alpha);
        grad1{1} = 2*lambda_scale*Ht{1}((H{1}(u1t{1},u1t{2},a) - Yt_alpha{1}),a);
        grad1{2} = 2*lambda_scale*Ht{2}((H{2}(u1t{1},u1t{2},a) - Yt_alpha{2}),a);
        grad1{3} = 2*lambda_scale*Ht{3}((H{3}(u1t{3},u1t{4},a) - Yt_alpha{3}),a);
        grad1{4} = 2*lambda_scale*Ht{4}((H{4}(u1t{3},u1t{4},a) - Yt_alpha{4}),a);
        
        gu{1} = (grad1{1}(1,:) + grad1{2}(1,:))./2;
        gu{2} = (grad1{1}(2,:) + grad1{2}(2,:))./2;
        gu{3} = (grad1{3}(1,:) + grad1{4}(1,:))./2;
        gu{4} = (grad1{3}(2,:) + grad1{4}(2,:))./2;

for i = 1:4
        grad2{i} = nuo*(u1t{i} - u2t{i}); % + mu*(u1t{i} - phi{i});
        grad{i} = gu{i}.' + grad2{i};
        g{i} = u1t{i} - gamma1(i)*grad{i};
        % proximity step
        vr{i} = min(max(real(g{i}), theta_minoR(:,i)), theta_maxoR(:,i));
        vi{i} = min(max(imag(g{i}), theta_minoI(:,i)), theta_maxoI(:,i));
        u1t{i} = vr{i} + 1i*vi{i};
end
    end
    % ----------------------------------------------------
    for i = 1:4
%         Ut1__ = u1t{i};
        Ut1__(i,:,a) =  flipud(u1t{i});
    end
    id(a) = true;
end

end

