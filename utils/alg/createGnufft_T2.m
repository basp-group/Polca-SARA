function [G, ll, v, V] = createGnufft_T2(D1, D2, om, K, S, J, W)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% square images and Fourier planes, to be modified
%
% om : {T,1} contain the frequencies to center the convolved kernels
% U1 : {T,1}, U1{t} [na_t|S^2]
% -> to be applied to each time instant independently
%-------------------------------------------------------------------------%
% Remarks:
%
% - keep frequency centers in memory [M|1] ?
% - keep supports in memory ? [M|2*S^2-1]
% - compare conv2 and doubleloop solution for small kernels -> conv2 quicker.
%-------------------------------------------------------------------------%
% Code: P.-A. Thouvenin, [08/12/2017]
% [Debug: check indices and values from the 1D NUIFFT] 
%-------------------------------------------------------------------------%
%%
Q  = 2*S - 1;             % size of the kernels after convolution
Q2 = Q^2;
J2 = J^2;

Qprime = floor(Q/2);      % Q is always odd (one possibility only)
tmp1 = (-Qprime:Qprime).';
tmp1 = tmp1(:, ones(Q,1));

[~, na] = size(D1{1}); % [S2, na, T] number of antennas at time t
M_true  = na*(na-1)/2; % number of acquisitions -> check value...
M = size(om, 1);    % M = M_true if all the measurements are present, M < M_true if some of the data have been flagged
% v = zeros(Q^2,M_true);   % convolution values (stored in column for each pair)
m = 0;              % global counter

%% Perform 2D convolutions and gridding using D1 and D2 (to be performed possibly in parallel)
for alpha = 1:na-1
    for beta = alpha+1:na % modify the double loop to exclusively select the appropriate elements, apply nonzeros on W
        % 2D convolutions
         m =m+1;
         q= 0;
        for i =1:4
            for j =1:4
                q = q+1;
%         v(:,q) = reshape(conv2(rot90(reshape(D1(:,alpha),[S,S]),2),reshape(conj(D2(:,beta)),[S,S])), [Q^2,1]); % only select the appropriate entries...
v{q}(:,m) = reshape(conv2(rot90(reshape(D1{i}(:,alpha),[S,S]),2),reshape(conj(D2{j}(:,beta)),[S,S])), [Q^2,1]); % only select the appropriate entries...
            end
        end
    end 
end

% Generate indices in the sparse G matrix
if rem(J,2) > 0 %odd
   c0 = round(om.*K/(2*pi)) - (J+1)/2; % [M, 2]
else
   c0 = floor(om.*K/(2*pi)) - J/2; 
end
kdy = bsxfun(@plus, (1:J).', c0(:,1).'); % [J M]
kdx = bsxfun(@plus, (1:J).', c0(:,2).'); % [J M]
ii = mod(bsxfun(@plus, tmp1(:), reshape(kdy, [1,J,M])), K(1)) + 1;  % [Q2,J,M] %row indices of the elements within each area, 
                                                                   % whose leftmost element row indices are given above
jj = mod(bsxfun(@plus, reshape(tmp1.', [Q2,1]), reshape(kdx, [1,J,M])), K(2)) + 1; % [Q2, J, M] % column indices ...
ll_ = reshape(bsxfun(@plus, reshape((jj-1)*K(1), [Q2,1,J,M]), reshape(ii, [Q2,J,1,M])), [Q2,J2,M]);
ll = repmat(ll_,[1,4,1]); %[Q2,4*J2,M]

% Duplicate values to have all the convolutions centered in the different elements
for qq=1:q
    V_ = bsxfun(@times, reshape(v{qq}, [Q2, 1, M_true]), reshape(W, [1, J2, M_true])); % [Q2, J2, M] 
    V_((isnan(V_(:)))) = [];
    V(qq,:,:,:) = V_;
end
% V = bsxfun(@times, reshape(v, [Q2, 1, M_true]), reshape(W, [1, J2, M_true])); % [Q2, J2, M]   

V1 = cat(3,V(1,:,:,:),V(2,:,:,:),V(9,:,:,:),V(10,:,:,:));
V2 = cat(3,V(3,:,:,:),V(4,:,:,:),V(11,:,:,:),V(12,:,:,:));
V3 = cat(3,V(5,:,:,:),V(6,:,:,:),V(13,:,:,:),V(14,:,:,:));
V4 = cat(3,V(7,:,:,:),V(8,:,:,:),V(15,:,:,:),V(16,:,:,:));


% V(isnan(V(:))) = []; % [J2*M, 1] there are zeros in W at the positions where the 
                 % measurements are missing, right size once the zeros are 
                 % filtered
% V1 = V{1}(     
% Generate row indices (within G) [remark: Jt^2 nz elements per row]
kk = bsxfun(@times, ones(J2*Q2,1), 1:4*M); % [J2*Q2,4*M]
% try


G{1} = sparse(kk(:), ll(:), V1(:), 4*M, K(1)*K(2));  % [M|K(1)*K(2)]
G{2} = sparse(kk(:), ll(:), V2(:), 4*M, K(1)*K(2));  % [M|K(1)*K(2)]
G{3} = sparse(kk(:), ll(:), V3(:), 4*M, K(1)*K(2));  % [M|K(1)*K(2)]
G{4} = sparse(kk(:), ll(:), V4(:), 4*M, K(1)*K(2));  % [M|K(1)*K(2)]

% catch 
%     keyboard
% end
end
