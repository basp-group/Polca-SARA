function Ua = computeUa2(Da, F, T, P, Gt, scale) % idt
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% [../../..]
%
%-------------------------------------------------------------------------%
%%
% [S2, T] = size(Da);
% [S2, P] = size(Ua);
% Compute D1 (or D2) for the 2D convolutions to be performed later
% Ua = so_ifft_adj((Da*conj(Gt)).', P, F, scale).'; % [S2, P] % zero entries are filtered from the transform 
%%- Ua = so_ifft_adj((Da(:, id_t)*conj(Gt(id_t,:))).', P, F, scale).';

% Ua = so_ifft_adj(Da.', T, P, scale).'; % [S2, P]

% Corrected version
% c = floor(T/2) + 1;
% p = floor(P/2);
% y = fftshift(fft(Da.', T, 1), 1); % [T, S2]
% Ua = y(c-p:c+p, :).'; % [S2, P]

% Corrected version
c = floor(F/2) + 1;
p = floor(P/2);
y = fftshift(fft(Da.', F, 1),1)/sqrt(F); % [T, S2]
Ua = y(c-p:c+p, :).'; % [S2, P]

end
