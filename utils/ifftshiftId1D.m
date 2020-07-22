function id_shift = ifftshiftId1D(Nfft,id)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

id1 = (id <= round(Nfft/2));
id_shift = zeros(size(id));
id_shift(id1) = id(id1) + floor(Nfft/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_shift(~id1) = id(~id1) - round(Nfft/2);

end

