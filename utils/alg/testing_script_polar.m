clc; clear all; close all;
format compact;
% [16/11/2017]
% [05/12/2017] Last debug: - the problem lied in the way I tested the results
% (need for symmetric x_hat and U1 = U2).
%                          - works currently when J and S are odd (see modifications to have it working with J or S even)
%                          - modify neighbourhood computation (consistency with Fessler's procedure and with PURIFY... to be checked)

% CONCLUSION: - keep the second version (create X_ab, no D necessary) [i.e., the shifts only affect the X matrix]
% - the results are not systematically onsistent with the first version,
% contrary to the second

% [08/12/2017] Add temporal prior (symmetry is preserved, the values coincide between the two formulations)
% -> check if the temporal NUFFT is correctly applied (value of nshift, ...)
% -> think about a parallel implementation of the solver

addpath utils
addpath src
addpath lib/nufft
addpath lib
addpath synth
rng(3);

%% Definition of the variables

im_choice = 'cyg_a'; % 'avery'; % W28_256 M31
cov_type = 'vlaa'; % vlaa, meerkat
test_number = '1';
test_aux = 'nu10';
S = 3;

rng(3);

% Loading data and initialization variables
name = [cov_type, '_', im_choice, '_', test_number,'_',num2str(S)];
 load(['synth_data_reg_eusipco_', name,'.mat'])

D1 = D_th;
D2 = D1;


% Test H1a, H2a, direct-operator
for i=1:4
    x_hat{i} = reshape(fftshift(reshape(so_fft2(b{i}, K, scale), K)), [prod(K), 1]); % FT of brightness matrix
end

H1a = zeros(8, na-1, S2, T);
H2a = zeros(8, na-1, S2, T);
Y11 = zeros(na, na, T);
Y22 = zeros(na, na, T);
id = true(na,1);
for a = 1:na
    id(a) = false;
    for t = 1:T
        % H1
        id_a = find(~isnan(squeeze(Y0(a, id, t)))); % find
        om1 = squeeze(Omega{t}(a,id,1)).'; % [S2^2*J2*(na-1), na]
        om1(isnan(om1)) = [];
        om2 = squeeze(Omega{t}(a,id,2)).'; % check dimension
        om2(isnan(om2)) = [];
        id1 = indices4Xhat(S, J, K, [om1, om2]);
        if ~isempty(id_a)
%             D2_at = computeD(U2(:, id, :), T, Gt(t, :), scale, T); % [S2, na-1, 1] -> problem with this option
            D2_at = D2(:, id, t);
            H1a(id_a, :, t) = createH1a(x_hat2(id1), squeeze(V2{t}(:, a, id)), D2_at, J, S, id_a).'; % [n, S2]    
        end
        Y1(a, id, t) = (flipud(D1(:,a,t)).'*H1a(:,:,t).'); % (id, a, t)
        
%         keyboard
        
        % H2
        id_a = find(~isnan(squeeze(Y0(id, a, t)))); % find
        om1 = squeeze(Omega{t}(id,a,1)); % [S2^2*J2*(na-1), na]
        om1(isnan(om1)) = []; % check format
        om2 = squeeze(Omega{t}(id,a,2));
        om2(isnan(om2)) = [];
        id1 = indices4Xhat(S, J, K, [om1, om2]);
        if ~isempty(id_a)
%             D1_at = computeD(U1(:, id, :), T, Gt(t, :), scale, T); % [S2, na-1, 1]
            D1_at = D1(:, id, t);
            H2a(id_a, :, t) = createH2a(x_hat2(id1), squeeze(V2{t}(:, id, a)), D1_at, J, S, id_a); % [n, S2]    
        end
        Y2(id, a, t) = (H2a(:,:,t)*conj(D2(:,a,t))); % (a, id, t)
    end
    
    % Test direct operator
    Y11(a, id , :) = reshape(direct_operator(flipud(squeeze(U1(:, a, :))), H1a, F, T, Gt, scale, []), [1, na-1, T]);
    Y22(id, a , :) = reshape(direct_operator(fliplr(conj(squeeze(U2(:, a, :)))), H2a, F, T, Gt, scale, []), [na-1, 1, T]);
    
    id(a) = true;
end
err1 = norm(Y1(:) - Y0(:))
err11 = norm(Y11(:) - Y0(:))
err2 = norm(Y2(:) - Y0(:))
err22 = norm(Y22(:) - Y0(:))



% Test with fliplr (formulation for updateU2)
Ha = ones(na-1, S*S, T);
x0 = randn(S*S, P);
Hx = direct_operator(fliplr(x0), Ha, F, T, [], ones(T, 1), []);
y0 = randn(na-1, T);
Hty = fliplr(adjoint_operator(y0, Ha, F, T, Gt, ones(T, 1), P));
a = sum(sum(Hty.*x0));
b = sum(sum(Hx.*y0)); % a = conj(b)
norm(b - conj(a))

% Test matrix multiplication
a = randn(na-1,S2,T);
b = 2*randn(S2, T);
c = zeros(na-1, T);
for t = 1:T
   c(:, t) = a(:, :, t)*b(:, t);
end
c2 = squeeze(sum(bsxfun(@times, a, reshape(b, [1, S2, T])), 2));
norm(c(:)-c2(:))

% Test compute D and compute Da
Da1 = computeDa(squeeze(U1(:, 1, :)), F, Gt, scale, T, []);
Da2 = squeeze(computeD(U1(:, 1, :), F, Gt, scale, T));
isequal(Da1, Da2)

% Test compute Da (test similarity etween computeD ad computeDa)
id(1)=  false;
D1_at = computeD(U1(:, id, :), F, [], scale, T); % [S2, na-1, 1]
D1_at_ref = D1(:, id, :);
isequal(D1_at, D1_at_ref)




%%

% Test with fliplr (formulation for updateU2)
Ha = ones(8,na-1, S*S, T);
u1 = randn(S*S, P,4);
u2 = u1;
Hx = direct_operator(fliplr(x0), Ha, F, T, [], ones(T, 1), []);
y0 = randn(na-1, T);
Hty = fliplr(adjoint_operator(y0, Ha, F, T, Gt, ones(T, 1), P));
a = sum(sum(Hty.*x0));
b = sum(sum(Hx.*y0)); % a = conj(b)
norm(b - conj(a))


Ht{1} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P,1);
Ht{3} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P, 1);
Ht{2} = @(y) adjoint_operator2(y, h11, h12, F, T, Gt, scale_t, P, 1);
Ht{4} = @(y) adjoint_operator2(y, h13, h14, F, T, Gt, scale_t, P, 1);

