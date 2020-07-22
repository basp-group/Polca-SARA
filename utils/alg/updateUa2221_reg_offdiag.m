function Ua2__ = updateUa2221_reg_offdiag(Ya, Ua1_, Ua2_, H2a, Gt, scale_t, F, T, JU2o, nuo, mu, lambda_scale, theta_minoR, theta_maxoR, theta_minoI, theta_maxoI)

% Update the DDE term 2 related to antenna a (PALM step).
%-------------------------------------------------------------------------%
% Input:
% > Ya           : column a of the data matrix Y with redundancy [na-1, T]
% > Ua2          : DDEs related to antenna a [S2, P]
% > Ua1          : DDEs related to antenna a [S2, P]
% > H2a          : linear operator involved in the update step [S2, na - 1]
% > Gt           : (unused) gridding matrix for temporal NUFFT
% > scale_t      : (unused) scaling coeff. for temporal NUFFT
% > F            : size of the temporal Fourier space
% > T            : number of time instants
% > JU2o         : maximum number of iterations
% > nuo          : regularization parameter (distance between U1a and U2a)
% > mu           : regularization parameter (distance between U1a and Phi)
% > lambda_scale : = 1/prod(K), K size of the spatial Fourier space
% > theta_minoR, theta_maxoR : box constraints on the real part of U1a
% > theta_minoI, theta_maxoI : box constraints on the imaginary part of U1a
%
% Output:
% < Ua2__          : DDEs related to antenna a [S2, P]
%
%-------------------------------------------------------------------------%
%% 
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%
%%

% U1: [S2, na, P] -> D1: [S2, na, T] 
% U2: [na, S2, P] -> D2: [na, S2, T] 
% Y: [na, na, T] -> Ya: [na-1, T]
% Ua1: [S2, P]


[S2, P,~] = size(Ua2_);
id_a = find(Ya{1} == 0); % indices of the missing measurements (represented by 0) 
phi = sparse(floor(S2/2)+1, floor(P/2)+1, 1, S2, P);
up_u2 = 1;

% Function handles to create H matrix for each measurement made by the
% correlator
% H{1} = @(u1,u2,a) (H2t{1,a}+ H2t{2,a})*u1 + (H2t{3,a}+ H2t{4,a})*u2;
% H{3} = @(u1,u2,a) (H2t{5,a}+ H2t{6,a})*u1 + (H2t{7,a}+ H2t{8,a})*u2;
% H{2} = H{1};
% H{4} = H{3};
% 
% 
z = size(H2a);
h11 = reshape((H2a(1,:,:,:) + H2a(2,:,:,:)), [z(2:end) 1]);
h12 = reshape((H2a(3,:,:,:) + H2a(4,:,:,:)), [z(2:end) 1]);
h13 = reshape((H2a(5,:,:,:) + H2a(6,:,:,:)), [z(2:end) 1]);
h14 = reshape((H2a(7,:,:,:) + H2a(8,:,:,:)), [z(2:end) 1]);

% h11 = squeeze(H2a(1,:,:,:) + H2a(2,:,:,:));
% h12 = squeeze(H2a(3,:,:,:) + H2a(4,:,:,:));
% h13 = squeeze(H2a(5,:,:,:) + H2a(6,:,:,:));
% h14 = squeeze(H2a(7,:,:,:) + H2a(8,:,:,:));

H{1} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h11, h12, F, T, Gt, scale_t, id_a, S2); % do not forget the fliplr (property F(x*)[n] = {F(x)*}[-n]) 
H{3} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h13, h14, F, T, Gt, scale_t, id_a, S2);
H{2} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h11, h12, F, T, Gt, scale_t, id_a, S2);
H{4} = @(u1,u2) direct_operator2(fliplr(u1),fliplr(u2), h13, h14, F, T, Gt, scale_t, id_a, S2);



% Adjoint operators
% Ht{1} = @(y,a) [(H2t{1,a}+H2t{2,a})'*y, (H2t{3,a}+H2t{4,a})'*y].';
% Ht{3} = @(y,a) [(H2t{5,a}+H2t{6,a})'*y, (H2t{7,a}+H2t{8,a})'*y].';
% Ht{2} = Ht{1};
% Ht{4} = Ht{3};
% 
% 
Ht{1} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P, S2, up_u2);
Ht{3} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P, S2, up_u2);
Ht{2} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P, S2, up_u2);
Ht{4} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P, S2, up_u2);


 
for i = 1:4
    Ua1{i} = squeeze(Ua1_(:,:,i)); %reshape(Ut1_(i,:,:),[z1(2:end) 1]);
    Ua2{i} = squeeze(Ua2_(:,:,i)); %reshape(Ut2_(i,:,:),[z2(2:end) 1]);
    u1a{i} = conj(Ua1{i}); % [S2, P]
    u2a{i} = conj(Ua2{i});
    Lips_temp(i) = nuo + 2*lambda_scale*pow_method_stokes_cal_reg(@(x1,x2) H{i}(x1,x2), @(x) Ht{i}(x), size(u1a{i})); % + muo % /F
    gamma2(i) = 1.9/Lips_temp(i) ;
end

% ----------------------------------------------------
% Iterations
% ----------------------------------------------------
    for q = 1:JU2o
        % gradient step
        %         grad1 = 2*lambda_scale*Ht2(:, :, a)'*(Ht2(:, :, a)*ut2 - Yt_alpha);
        grad1{1} = 2*lambda_scale*Ht{1}((H{1}(u2a{1},u2a{3}) - Ya{1})); %/F; %u2t{2}
        grad1{2} = 2*lambda_scale*Ht{2}((H{2}(u2a{2},u2a{4}) - Ya{2})); %/F; %u2t{3}
        grad1{3} = 2*lambda_scale*Ht{3}((H{3}(u2a{1},u2a{3}) - Ya{3})); %/F; %u2t{2}
%         grad1{4} = 2*lambda_scale*Ht{4}((H{4}(u2a{3},u2a{4}) - Ya{4})); %/F; 
%         
%         gu{1} = (grad1{1}(:,:,1) + grad1{3}(:,:,1))./2;
%         gu{2} = (grad1{1}(:,:,2) + grad1{3}(:,:,2))./2;
%         gu{3} = (grad1{2}(:,:,1) + grad1{4}(:,:,1))./2;
%         gu{4} = (grad1{2}(:,:,2) + grad1{4}(:,:,2))./2;
%         

%%%% Rectified
        grad1{4} = 2*lambda_scale*Ht{4}((H{4}(u2a{2},u2a{4}) - Ya{4})); %/F; 
        
        gu{1} = (grad1{1}(:,:,1) + grad1{3}(:,:,1))./2;
        gu{3} = (grad1{1}(:,:,2) + grad1{3}(:,:,2))./2;
        gu{2} = (grad1{2}(:,:,1) + grad1{4}(:,:,1))./2;
        gu{4} = (grad1{2}(:,:,2) + grad1{4}(:,:,2))./2;
        
        
        for i = 2:3 %1:4
            grad2{i} = nuo*(u2a{i} - u1a{i}); % + mu*(u2t{i} - phi{i});
            grad{i} = gu{i} + grad2{i};
            g{i} = u2a{i} - gamma2(i)*grad{i};
            % proximity step
            vr{i} = min(max(real(g{i}), theta_minoR(:,:,i)), theta_maxoR(:,:,i));
            vi{i} = min(max(imag(g{i}), theta_minoI(:,:,i)), theta_maxoI(:,:,i));
            u2a{i} = vr{i} + 1i*vi{i};
        end
        
    end
    % ----------------------------------------------------
    for i = 1:4
%         Ut2__ = u2t{i}; 
         Ua2__(:,:,i) = conj(u2a{i});
    end

end

