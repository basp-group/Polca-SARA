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
na = 5;
T = 10;
om_alpha = rand(na,2,T)*pi/10;
t_samples = cumsum(rand(T,1)*pi/20); % see time scale for the time samples..., see if no boundary effects

M = na*(na-1)/2;
N = [128, 128]; % mandatory choice: even value (K and N)
K = 2*N;
om = zeros(M, 2, T);
P = 5;    % mandatory : odd value
F = 2048; % 2048
Jt = 3;

for t = 1:T
    q = 0;
    for a = 1:na-1
        for b = a+1:na
            q = q+1;
            om(q,:,t) = om_alpha(a,:,t) - om_alpha(b,:,t);
        end
    end
end

J = 3;  % mandatory choice: odd value (modify code to allow even values? to be seen...)
S = 3;  % mandatory choice: odd value
S2 = S^2;
J2 = J^2;
Js = J + S - 1;
Js2 = Js^2;
Q = 2*S - 1;
Q2 = Q^2;

x = fitsread('M31.fits');
x = imresize(x, N);
x = (x+abs(x))./2;
x = x./max(x(:));

sig = 5e-2;
c = floor(S/2)*S + floor(S/2) + 1 ; % center of the support (local index) [S2, na, P]
p = floor(P/2) + 1;
U1 = (sig*randn(S2, na, P,4) + 1i*sig*randn(S2, na, P,4)); % with temporal prior
U1(c, :, p) = sqrt(F);
U2 = U1;

keyboard

% Duplicate entries for create H1 and H2
% Spatial gridding coefficients and associated frequencies (create Y at the same time ? see formats...)
V = cell(T, 1);
V2 = cell(T, 1);
Omega = cell(T, 1);
indices = getIndicesLow(na); % indices for the lower part of Y

for t = 1:T    
    
    % frequencies
    w_temp = zeros(na);
    w_temp(indices) = om(:,1,t); % to be seen (which variable should I use?)
    Omega{t}(:,:,1) = -w_temp + w_temp.';
    w_temp(indices) = om(:,2,t);
    Omega{t}(:,:,2) = -w_temp + w_temp.'; % [na, na, 2] for each t
    
    % spatial gridding coefficients
    st = compute_interp_coeffs(om(:,:,t), N, [J, J], K, N/2); % [Ny/2, Nx/2], N = K/2 [Mt, J2] 
    V{t} = st.uu; % [J2, M] 
    v_temp = flipud(conj(st.uu)); % revoir...
    W = zeros(na, na, J2);
    w_temp = zeros(na);
    for j = 1:J2
        w_temp(indices) = st.uu(j,:); % explained by the symmetry of the coefficients
        w_temp = w_temp.';
        w_temp(indices) = v_temp(j,:);
        W(:, :, j) = w_temp;
    end
    V2{t} = permute(W, [3, 1, 2]); % [J2, na, na]
    sp_scale_prev = st.sn;
%     if t > 1
%         isequal(st.sn, sp_scale_prev)
%     end
end

% Create D1 and D2
% [Gt, scale] = createGt(t_samples, Jt, P, F, T);
Gt = []; scale = [];
D1 = computeD(U1, F, Gt, scale, T); % [S2,na,T]
D2 = computeD(U2, F, Gt, scale, T); % [S2,na,T]

Y = zeros(na, na, T);
Y0 = zeros(na, na, T);
Y1 = zeros(na, na, T);
Y2 = zeros(na, na, T);
y = zeros(M, T);

% Data generation
for t = 1:T
    % Create operator G [for each time instant]
    G = createGnufft_T2(D1(:, :, t), D2(:, :, t), [om(:, 1, t), om(:, 2, t)], K, S, J, V{t});
    A = @(x) G*so_fft2(x, K, st.sn); % voir si st.sn varie en fonction de t
    y(:, t) = A(x);
%     y(:,t) = G*x_hat; % [M, 1]
    Y0(:,:,t) = createY(y(:,t), na);
end

% Test H1a, H2a, direct-operator
x_hat2 = reshape(fftshift(fft2(x.*st.sn, K(1), K(2))), [prod(K), 1]);
H1a = zeros(na-1, S2, T);
H2a = zeros(na-1, S2, T);
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

