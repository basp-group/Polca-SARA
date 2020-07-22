function  [H2a1,H2a2,H2a3,H2a4,H2a5,H2a6,H2a7,H2a8] = createH2a_32_reg(X1, X2, X3, X4, V, U1, J, S, id_nnz)
% Compute the linear operator Ha involved in the update of the second
% DDE term for a single antenna a (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > X     : part of xhat [S2*S2*J2*(na - 1), 1]
% > V     : part of the spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na - 1]
% > D1_at : D1(:, :, t), where the column a has been removed 
%           D1 : [S2, na, T] -> D1_at : [S2, na-1]
% > J     : size of the gridding kernels
% > S     : size of the DDE support (spatial Fourier domain)
% > id_nnz: position of missing measurements (temporal dimension)
%
% Output:
% < H2a   : linear operator Ha [S2, na - 1]
%
%-------------------------------------------------------------------------%


S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 
z = size(U1);

for i =1:4
D1_at{i} = reshape(U1(i,:,:,:),[z(2:end) 1]);  %squeeze(U1(i,:,:));
end

for q = 1:n
    Xq{1} = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X1((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
    Xq{2} = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X2((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
    Xq{3} = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X3((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
    Xq{4} = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X4((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified

    H2a1(q,:) = flipud(D1_at{1}(:,q)).'*Xq{1}; % select the appropriate portion of D2_at to account for missing measurements
    H2a2(q,:) = flipud(D1_at{2}(:,q)).'*Xq{3}; 
    H2a3(q,:) = flipud(D1_at{1}(:,q)).'*Xq{2}; 
    H2a4(q,:) = flipud(D1_at{2}(:,q)).'*Xq{4}; 
    H2a5(q,:) = flipud(D1_at{3}(:,q)).'*Xq{1}; 
    H2a6(q,:) = flipud(D1_at{4}(:,q)).'*Xq{3}; 
    H2a7(q,:) = flipud(D1_at{3}(:,q)).'*Xq{2}; 
    H2a8(q,:) = flipud(D1_at{4}(:,q)).'*Xq{4}; 
end

% 
 end
% 
