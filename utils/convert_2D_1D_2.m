function n = convert_2D_1D_2(c, K)

% c : position 2D dans {-K1/2, ..., K1/2-1}x{-K2/2, ..., K2/2-1}
% n : position 1D dans {1,...,K1*K2}
% N : taille de l image = K1xK2

n = (c(:,2) + floor(K(2)/2))*K(1) + c(:,1) + floor(K(1)/2) +1;

end