function Il = getIndicesLow(N)
%getIndicesLow
%   Get the linear indices corresponding to the lower part of a square
%   matrix of size N (construction of Y)

Il = zeros(N*(N-1)/2,1);
k = 1;
for j = 1:N-1
    Il(k:k+N-j-1) = (j-1)*N + (j+1:N)';
    k = k + N - j;
end

% Iu = N^2 - flipud(Il) + 1; % sorted indices coreesponding to the upper part of the matrix

end

