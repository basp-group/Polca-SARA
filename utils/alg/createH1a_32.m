function [H1a1,H1a2,H1a3,H1a4,H1a5,H1a6,H1a7,H1a8] = createH1a_32(x_hat, V, U2, J, S, id_nnz,id1,id, zz2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% om_at : same construction as Y
% check similar construction for V (avoid flip and conj?)
% D2 : [S2, na, T] -> D2_at : [S2, na-1]
% om_at : [na-1, 2]

S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 

% H1a = zeros(S2, n);
% for q = 1:n
%     Xq = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
% %     Xq = sum(bsxfun(@times, reshape(V(:,q), [1, 1, J2]), reshape(X{q}, [S2, S2, J2])), 3); % to be verified, currently error with reshape
%     H1a(:,q) = Xq*conj(D2_at(:,q)); % select the appropriate portion of D2_at to account for missing measurements
%     % H1a(:,q) = Xq*conj(D2_at(:,id_nnz(q)));
% end

% keyboard

% V: [J2, n]
% Xq = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
% H1a = reshape(sum(bsxfun(@times, permute(Xq, [2, 1, 3,]), reshape(conj(D2_at), [S2, 1, n])), 1), [S2, n]); % [S2, n]

% X = x_hat{1}(id1);
% D2_at = reshape(U2{t}(1,:,id),[zz2{t}(2:end) 1]);

for i =1:4
    X = x_hat{i}(id1);
    Xq{i} = reshape(sum(bsxfun(@times, reshape(V, [J2, 1, 1, n]), permute(reshape(X, [S2, S2, J2, n]), [3, 1, 2, 4])), 1), [S2, S2, n]); % [S2, S2, n]
    D2_at{i} = reshape(U2(i,:,id),[zz2(2:end) 1]);
end


H1a1 = (reshape(sum(bsxfun(@times, permute(Xq{1}, [2, 1, 3]), reshape(conj(D2_at{1}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]'
H1a2 = (reshape(sum(bsxfun(@times, permute(Xq{3}, [2, 1, 3]), reshape(conj(D2_at{1}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]
H1a3 = (reshape(sum(bsxfun(@times, permute(Xq{2}, [2, 1, 3]), reshape(conj(D2_at{3}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{2}
H1a4 = (reshape(sum(bsxfun(@times, permute(Xq{4}, [2, 1, 3]), reshape(conj(D2_at{3}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{2}
H1a5 = (reshape(sum(bsxfun(@times, permute(Xq{1}, [2, 1, 3]), reshape(conj(D2_at{2}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{3}
H1a6 = (reshape(sum(bsxfun(@times, permute(Xq{3}, [2, 1, 3]), reshape(conj(D2_at{2}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n] %D2_at{3}
H1a7 = (reshape(sum(bsxfun(@times, permute(Xq{2}, [2, 1, 3]), reshape(conj(D2_at{4}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]
H1a8 = (reshape(sum(bsxfun(@times, permute(Xq{4}, [2, 1, 3]), reshape(conj(D2_at{4}), [S2, 1, n])), 1), [S2, n])).'; % [S2, n]


end
