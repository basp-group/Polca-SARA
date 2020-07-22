function y = direct_operator2(x1,x2,Ha1,Ha2,F,T,Gt,scale,id_a, S2)

% x1: [S2,P]
% x2: [S2,P]
% Ha1: [ na-1,  T]
% Ha2: [ na-1, T]

y1 = squeeze(sum(bsxfun(@times, Ha1, reshape(computeDa2(x1, F, Gt, scale, T, id_a), [1, S2, T])), 2)); % [na-1, T2]
y2 = squeeze(sum(bsxfun(@times, Ha2, reshape(computeDa2(x2, F, Gt, scale, T, id_a), [1, S2, T])), 2)); % [na-1, T2]

y = y1 + y2;
end