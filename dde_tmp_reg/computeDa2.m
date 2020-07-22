function Da = computeDa2(Ua, F, Gt, scale, T, id_a)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% [../../..]
%
%-------------------------------------------------------------------------%
%%
% Version with Temporal NUIFFT
% D: [T, na*S2]
% U: [P, na*S2]
% [S2, P] = size(Ua);
% Compute D1 (or D2) for the 2D convolutions to be performed later
% Da = (Gt(id_t,:)*so_ifft(Ua.', F, scale)).'; % [S2, T'] % zero entries are filtered from the transform
% Da = (Gt*so_ifft(Ua.', F, scale)).'; % [S2, T]
% Da = so_ifft(Ua.', T, F, scale);
% Da = Da(1:T, :).'; % [S2, T]

% % Corrected version
% [S2, P] = size(Ua);
% x = zeros(T, S2);
% c = floor(T/2) + 1;
% p = floor(P/2);
% x(c-p:c+p, :) = Ua.';
% Da = T*ifft(ifftshift(x, 1), T, 1).'; % [S2, T]

% Test with additional 0-padding
[S2, P] = size(Ua);
x = zeros(F, S2);
c = floor(F/2) + 1;
p = floor(P/2);
x(c-p:c+p, :) = Ua.';
Da = sqrt(F)*ifft(ifftshift(x, 1), F, 1).'; % [S2, T]
Da = Da(:, 1:T);
Da(id_a) = 0; % set to 0 the time instants corresponding to missing measurements

end
