function [I,Ic] = getSupportRect(K, S, N1, N2)
% Generate the indices corresponding to the pixels contained in the 
% frequency region defined by the support S, centered in the spatial
% frequencies given in K.
%
% Input
% > K : spatial frequencies, assumed ordered by t (Na contiguous 
%       elements correspond to a single time instant) [Na|1]
%       K(q,:) \in {-N/2, ..., N/2 - 1}^2 if N even,
%                  {-floor(N/2), ..., floor(N/2)}^2 if N odd
% > S : support size (assume square support)
% > N : image size (rectangle images supported)
% > Na : number of antenna.
%
% Output
% < I  : linear indices of the pixel contained in the frequency region of 
%        size S, centered in the elements of K (linear indices) [S^2*Na, T]
% < Ic : indices corresponding to the entries on I [S^2, 2, NaT]
%
% Index of the pixels belonging to the neighbourhood of the pixel specified
% in K(q,:) : I(~isnan(I(:,t)),t), q = (t-1)Na + m
%-------------------------------------------------------------------------%
%%
% Remark:
% entries to be filtered out with ~isnan(), output format to be modified;
%-------------------------------------------------------------------------%
%%
% Code : Pierre-Antoine Thouvenin, October 9th 2017.
%%
% Auxiliary variables (to handle the cases where S is even or odd)
r = rem(S,2) - 1;
S_prime = floor(S/2);
Na = size(K,1);

% Center the elements in [1,N]^2
% K_prime = K + floor(N/2) + 1; % assuming the frequencies are already
% centered as appropriate

% Generate all the couples (3D)
tmp = (-S_prime:(S_prime - r))'; % row indices defining the elements in S
tmp2 = tmp(:,ones(S,1));         % column idices defining the elements in S
n_c = [reshape(tmp2, [S^2,1]), reshape(tmp2', [S^2,1])]; 
Ic = bsxfun(@plus, n_c, reshape(K', [1, 2, size(K, 1)])); % [S^2,2,Na]
Ic = reshape(permute(Ic, [1,3,2]), [S^2*Na,2]);           % [S^2*Na,2]
% keyboard

% Filter out invalid positions (i.e., out of the scope of the image)
id1 = bsxfun(@or, Ic(:,1) < 1, Ic(:,1) > N1);
id2 = bsxfun(@or, Ic(:,2) < 1, Ic(:,2) > N2);
Ic(id1,:) = NaN; % [S^2*Na, 2]
Ic(id2,:) = NaN; % [S^2*Na, 2]

% Convert to linear indices
I = reshape((Ic(:,2) - 1)*N2 + Ic(:,1), [S^2,Na]); % [S^2,Na]

% I(isnan(I)) = 0; % to be preserved by min/max operations
% Ic(id,:) = 0;
% renvoyer également les couples d'indices pour la création de matrice
% sparses ?
end
