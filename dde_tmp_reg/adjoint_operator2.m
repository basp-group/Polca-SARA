function x = adjoint_operator2(y, Ha1, Ha2, F, T, Gt, scale, P, S2, up_u2)

[n, ~] = size(y); % n = na - 1
x(:,:,1) = computeUa2(reshape(sum(bsxfun(@times, permute(conj(Ha1), [2, 1, 3]), ...
    reshape(y, [1, n, T])), 2), [S2, T]), F, T, P, Gt, scale); % [S2, P]

x(:,:,2) = computeUa2(reshape(sum(bsxfun(@times, permute(conj(Ha2), [2, 1, 3]), ...
    reshape(y, [1, n, T])), 2), [S2, T]), F, T, P, Gt, scale); % [S2, P]

if up_u2 == 1
    x(:,:,1) = fliplr(x(:,:,1));
    x(:,:,2) = fliplr(x(:,:,2));
end
end
