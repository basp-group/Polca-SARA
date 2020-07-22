function Ua1__ = updateUt122_reg(Ya, Ua1_, Ua2_, H1a, na_t, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here

% U2t: [S2, na]
% U1t: [S2, na]
% Yt: [na, na]
% H1t: [na-1, S2, na]

[S2, P, ~] = size(Ua1);
id_a = find(Ya{1} == 0); % indices of the missing measurements (represented by 0) 
phi = sparse(floor(S2/2)+1, floor(P/2)+1, 1, S2, P);

% Function handles to create H matrix for each measurement made by the
% % correlator
% H{1} = @(u1,u2,a) (H1t{1,a}+ H1t{3,a})*u1 + (H1t{2,a}+ H1t{4,a})*u2;
% H{2} = @(u1,u2,a) (H1t{5,a}+ H1t{7,a})*u1 + (H1t{6,a}+ H1t{8,a})*u2;
% H{3} = H{1};
% H{4} = H{2};

h11 = H1a(1,:,:,:) + H1a(3,:,:,:);
h12 = H1a(2,:,:,:) + H1a(4,:,:,:);
h13 = H1a(5,:,:,:) + H1a(7,:,:,:);
h14 = H1a(6,:,:,:) + H1a(8,:,:,:);

H{1} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a);
H{2} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a);
H{3} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a);
H{4} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a);

% Adjoint operators
% Ht{1} = @(y,a) [(H1t{1,a}+H1t{3,a})'*y, (H1t{2,a}+H1t{4,a})'*y].';
% Ht{2} = @(y,a) [(H1t{5,a}+H1t{7,a})'*y, (H1t{6,a}+H1t{8,a})'*y].';
% Ht{3} = Ht{1};
% Ht{4} = Ht{2};

Ht{1} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P);
Ht{2} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P);
Ht{3} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P);
Ht{4} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P);


for i = 1:4
    Ua1{i} = squeeze(Ua1_(:,:,i)); %reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ua2{i} = squeeze(Ua2_(:,:,i)); %reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1a{i} = flipud(Ua1{i});
    u2a{i} = flipud(Ua2{i});
    Lips_temp(i) = nuo + 2*lambda_scale*pow_method_stokes_cal(@(x1,x2) H{i}(x1,x2), @(x) Ht{i}(x), size(u1a{i})); % + muo
    gamma1(i) = 1.9/Lips_temp(i) ;
end

    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU1o
        % gradient step
%         grad1 = 2*lambda_scale*H1t(:, :, a)'*(H1t(:, :, a)*u1t - Yt_alpha);
        grad1{1} = 2*lambda_scale*Ht{1}(H{1}(u1a{1},u1a{2}) - Ya{1});
        grad1{2} = 2*lambda_scale*Ht{2}(H{2}(u1a{1},u1a{2}) - Ya{2});
        grad1{3} = 2*lambda_scale*Ht{3}(H{3}(u1a{3},u1a{4}) - Ya{3});
        grad1{4} = 2*lambda_scale*Ht{4}(H{4}(u1a{3},u1a{4}) - Ya{4});
        
        gu{1} = (grad1{1}(1,:) + grad1{2}(1,:))./2;
        gu{2} = (grad1{1}(2,:) + grad1{2}(2,:))./2;
        gu{3} = (grad1{3}(1,:) + grad1{4}(1,:))./2;
        gu{4} = (grad1{3}(2,:) + grad1{4}(2,:))./2;

for i = 1:4
        grad2{i} = nuo*(u1a{i} - u2a{i}); % + mu*(u1t{i} - phi{i});
        grad{i} = gu{i}.' + grad2{i};
        g{i} = u1a{i} - gamma1(i)*grad{i};
        % proximity step
        vr{i} = min(max(real(g{i}), theta_minoR(:,:,i)), theta_maxoR(:,:,i));
        vi{i} = min(max(imag(g{i}), theta_minoI(:,:,i)), theta_maxoI(:,:,i));
        u1a{i} = vr{i} + 1i*vi{i};
end
    end
    % ----------------------------------------------------
    for i = 1:4
%         Ut1__ = u1t{i};
        Ua1__(:,:,i) =  flipud(u1a{i});
    end

end

