function [D,Dt] = computeD_reg(U, F, Gt, scale, T)
% Compute the DDE kernels in the temporal domain.
%-------------------------------------------------------------------------%
% Input:
% > U     : DDEs in the Fourier domain (spatial and temporal) [S2, na, P, corr]
% > F     : size of the temporal Fourier domain
% > Gt    : temporal gridding matrix (unused, for possible extension to temporal NUFFT)
% > scale : scaling factor (unused, for possible extension to temporal NUFFT)
% > T     : number of time instants
%
% Output:
% < D     : DDEs in the temporal domain for each correlator [S2, na, T]
% < Dt    : DDEs in the temporal domain for each time instant [corr, S2, na]
%-------------------------------------------------------------------------%
%%
[S2, n, P,corr] = size(U);

x = zeros(F, S2*n);
c = floor(F/2) + 1;
p = floor(P/2);

for i =1:corr
%     UU = squeeze(U(:,:,:,i));
    UU = U(:,:,:,i);
    x(c-p:c+p, :) = reshape(permute(UU, [3, 1, 2]), [P, S2*n]); % [T, S2*n]
    D{i} = reshape(sqrt(F)*ifft(ifftshift(x, 1), F, 1).', [S2, n, F]); % T
    dd(i,:,:,:) = D{i}; % D{i}(:,:,1:T)
end

% To store values per time instant in a cell
Dt = cell(T,1);
parfor t =1:T
    Dt{t} = dd(:,:,:,t);
end

end
