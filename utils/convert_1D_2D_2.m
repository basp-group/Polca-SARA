function c = convert_1D_2D_2(n,K)

% n : position 1D dans {1,...,K1*K2}
% c : position 2D dans {-K1/2, ..., K1/2-1}x{-K2/2, ..., K2/2-1}
% N : taille de l image = K1xK2

d = floor(n/K(1)) ;
c = [n - (d*K(1)) - floor(K(1)/2) - 1, d - floor(K(2)/2)] ;

end