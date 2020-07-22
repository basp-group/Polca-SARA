function [H1a1,H1a2,H1a3,H1a4,H1a5,H1a6,H1a7,H1a8] = createH1a_32_reg(X1,X2,X3,X4, V, U2, J, S, id_nnz)
% Compute the linear operator Ha involved in the update of the first
% DDE term for a single antenna a (na antennas in total).
%-------------------------------------------------------------------------%
% Input:
% > X     : part of xhat [S2*S2*J2*(na - 1), 1]
% > V     : part of the spatial NUFFT gridding kernels with redundancy (similar to Y) [J2, na - 1]
% > D2_at : D2(:, :, t), where the column a has been removed 
%           D2 : [S2, na, T] -> D2_at : [S2, na-1]
% > J     : size of the gridding kernels
% > S     : size of the DDE support (spatial Fourier domain)
% > id_nnz: position of missing measurements (temporal dimension)
%
% Output:
% < H1a   : linear operator Ha [S2, na - 1]
%
%-------------------------------------------------------------------------%
% [15/03/2018], P.-A. Thouvenin.
%-------------------------------------------------------------------------%

S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 
z = size(U2);

Xq{1} = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X1, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
D2_at{1} = reshape(U2(1,:,:,:),[z(2:end) 1]);  % squeeze(U2(1,:,:));
Xq{2} = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X2, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
D2_at{2} = reshape(U2(2,:,:,:),[z(2:end) 1]);   %squeeze(U2(2,:,:));
Xq{3} = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X3, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
D2_at{3} = reshape(U2(3,:,:,:),[z(2:end) 1]);  % squeeze(U2(3,:,:));
Xq{4} = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X4, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
D2_at{4} = reshape(U2(4,:,:,:),[z(2:end) 1]);  % squeeze(U2(4,:,:));



H1a1 = (reshape(sum(bsxfun(@times, permute(Xq{1}, [2, 1, 3]), reshape(conj(D2_at{1}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]'
H1a2 = (reshape(sum(bsxfun(@times, permute(Xq{3}, [2, 1, 3]), reshape(conj(D2_at{1}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]
H1a3 = (reshape(sum(bsxfun(@times, permute(Xq{2}, [2, 1, 3]), reshape(conj(D2_at{3}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{2}
H1a4 = (reshape(sum(bsxfun(@times, permute(Xq{4}, [2, 1, 3]), reshape(conj(D2_at{3}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{2}
H1a5 = (reshape(sum(bsxfun(@times, permute(Xq{1}, [2, 1, 3]), reshape(conj(D2_at{2}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{3}
H1a6 = (reshape(sum(bsxfun(@times, permute(Xq{3}, [2, 1, 3]), reshape(conj(D2_at{2}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{3}
H1a7 = (reshape(sum(bsxfun(@times, permute(Xq{2}, [2, 1, 3]), reshape(conj(D2_at{4}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]
H1a8 = (reshape(sum(bsxfun(@times, permute(Xq{4}, [2, 1, 3]), reshape(conj(D2_at{4}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]


end
