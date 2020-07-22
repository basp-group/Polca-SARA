function H2a = createH2a_32(X, V, D1_at, J, S, id_nnz)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

% om_at : same construction as Y
% check similar construnction for V... (avoid flip and conj?)
% D1 : [S2, na, T] -> D2_at : [S2, na-1]
% om_at : [na-1, 2]

S2 = S^2;
J2 = J^2;
n = length(id_nnz); % n = na - 1 if all the measurements are present, n < na - 1 otherwise. 

H2a = zeros(n, S2);
for q = 1:n
    Xq = sum(bsxfun(@times, reshape(V(:,id_nnz(q)), [1, 1, J2]), reshape(X((q-1)*S2^2*J2+1:q*S2^2*J2), [S2, S2, J2])), 3); % to be verified
    H2a(q,:) = flipud(D1_at(:,q)).'*Xq; % select the appropriate portion of D2_at to account for missing measurements
end

end

