function Iu = getIndicesUp(N)
%getIndicesLow
%   Get the linear indices corresponding to the upper part of a square
%   matrix of size N.

Iu = zeros(N*(N-1)/2,1);
k = 1;
for j = 2:N
    Iu(k:k+j-2) = (j-1)*N + (1:j-1)';
    k = k + j - 1;
end

end

