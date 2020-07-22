function Ifft = getFourierSupport3(Ic, Nfft1, Nfft2, N1, N2)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here 

% to be merged with getSupport ?

% Convert to linear indices (w.r.t the image of size Nfft starting from the image of size N)
Ifft = zeros(size(Ic));
Ifft(:,1) = Ic(:,1) + floor(Nfft1/2) - floor(N1/2); % assume symmetric 0-padding around the image
Ifft(:,2) = Ic(:,2) + floor(Nfft2/2) - floor(N2/2);

end
