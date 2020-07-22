function [Ua1__] = updateUa1221_reg_offdiag(Ya, Ua1_, Ua2_, H1a, Gt, scale_t, F, T, JU1o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)
% Update the DDE term 1 related to antenna a (PALM step).
%-------------------------------------------------------------------------%
% Input:
% > Ya           : column a of the data matrix Y with redundancy [na-1, T]
% > Ua1          : DDEs related to antenna a [S2, P]
% > Ua2          : DDEs related to antenna a [S2, P]
% > H1a          : linear operator involved in the update step [S2, na - 1]
% > Gt           : (unused) gridding matrix for temporal NUFFT
% > scale_t      : (unused) scaling coeff. for temporal NUFFT
% > F            : size of the temporal Fourier space
% > T            : number of time instants
% > JU1o         : maximum number of iterations
% > nuo          : regularization parameter (distance between U1a and U2a)
% > mu           : regularization parameter (distance between U1a and Phi)
% > lambda_scale : = 1/prod(K), K size of the spatial Fourier space
% > theta_minoR, theta_maxoR : box constraints on the real part of U1a
% > theta_minoI, theta_maxoI : box constraints on the imaginary part of U1a
%
% Output:
% < Ua1__          : DDEs related to antenna a [S2, P]
%
%-------------------------------------------------------------------------%

% U1: [S2, na, P] -> D1: [S2, na, T] 
% U2: [S2, na, P] -> D2: [S2, na, T] 
% Y: [na, na, T] -> Ya: [na-1, T], with 0 for missing measurements
% Ua1: [S2, P]
% H1a: [na-1, S2, T]


[S2, P, ~] = size(Ua1_);
id_a = find(Ya{1} == 0); % indices of the missing measurements (represented by 0) 
phi = sparse(floor(S2/2)+1, floor(P/2)+1, 1, S2, P);

% Function handles to create H matrix for each measurement made by the
% % correlator
% H{1} = @(u1,u2,a) (H1t{1,a}+ H1t{3,a})*u1 + (H1t{2,a}+ H1t{4,a})*u2;
% H{2} = @(u1,u2,a) (H1t{5,a}+ H1t{7,a})*u1 + (H1t{6,a}+ H1t{8,a})*u2;
% H{3} = H{1};
% H{4} = H{2};

z = size(H1a);
h11 = reshape((H1a(1,:,:,:) + H1a(3,:,:,:)), [z(2:end) 1]);
h12 = reshape((H1a(2,:,:,:) + H1a(4,:,:,:)), [z(2:end) 1]);
h13 = reshape((H1a(5,:,:,:) + H1a(7,:,:,:)), [z(2:end) 1]);
h14 = reshape((H1a(6,:,:,:) + H1a(8,:,:,:)), [z(2:end) 1]);

% h12 = (H1a(2,:,:,:) + H1a(4,:,:,:));
% h13 = (H1a(5,:,:,:) + H1a(7,:,:,:));
% h14 = (H1a(6,:,:,:) + H1a(8,:,:,:));

H{1} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a, S2);
H{2} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a, S2);
H{3} = @(u1,u2) direct_operator2(u1,u2, h11, h12, F, T, Gt, scale_t, id_a, S2);
H{4} = @(u1,u2) direct_operator2(u1,u2, h13, h14, F, T, Gt, scale_t, id_a, S2);

% Adjoint operators
% Ht{1} = @(y,a) [(H1t{1,a}+H1t{3,a})'*y, (H1t{2,a}+H1t{4,a})'*y].';
% Ht{2} = @(y,a) [(H1t{5,a}+H1t{7,a})'*y, (H1t{6,a}+H1t{8,a})'*y].';
% Ht{3} = Ht{1};
% Ht{4} = Ht{2};

Ht{1} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P, S2, 0);
Ht{2} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P,S2,  0);
Ht{3} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P, S2,  0);
Ht{4} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P, S2,0);


for i = 1:4
    Ua1{i} = squeeze(Ua1_(:,:,i)); %reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ua2{i} = squeeze(Ua2_(:,:,i)); %reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1a{i} = flipud(Ua1{i});
    u2a{i} = flipud(Ua2{i});
    Lips_temp(i) = nuo + 2*lambda_scale*pow_method_stokes_cal_reg(@(x1,x2) H{i}(x1,x2), @(x) Ht{i}(x), size(u1a{i})); % + muo % /F
    gamma1(i) = 1.9/Lips_temp(i) ;
end

    % ----------------------------------------------------
    % Iterations
    % ----------------------------------------------------
    for q = 1:JU1o
        % gradient step
%         grad1 = 2*lambda_scale*H1t(:, :, a)'*(H1t(:, :, a)*u1t - Yt_alpha);
        grad1{1} = 2*lambda_scale*Ht{1}(H{1}(u1a{1},u1a{2}) - Ya{1}); %/F;
%         a1 = (H{4}(u1a{3},u1a{4})); 
        grad1{2} = 2*lambda_scale*Ht{2}(H{2}(u1a{1},u1a{2}) - Ya{2}); %/F;
        grad1{3} = 2*lambda_scale*Ht{3}(H{3}(u1a{3},u1a{4}) - Ya{3}); %/F;
        grad1{4} = 2*lambda_scale*Ht{4}(H{4}(u1a{3},u1a{4}) - Ya{4}); %/F;
        
        gu{1} = (grad1{1}(:,:,1) + grad1{2}(:,:,1))./2;
        gu{2} = (grad1{1}(:,:,2) + grad1{2}(:,:,2))./2;
        gu{3} = (grad1{3}(:,:,1) + grad1{4}(:,:,1))./2;
        gu{4} = (grad1{3}(:,:,2) + grad1{4}(:,:,2))./2;

for i = 2:3 %1:4
        grad2{i} = nuo*(u1a{i} - u2a{i}); % + mu*(u1t{i} - phi{i});
        grad{i} = gu{i} + grad2{i};
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

