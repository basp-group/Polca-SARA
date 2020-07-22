function id_shift = fftshiftId2D(I,J,id)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
% [26/10/17: ok]

% keyboard

id_2D_shift = zeros(numel(id),2); 
[idi,idj] = ind2sub([I,J], id(:));

idi1 = (idi <= floor(I/2));
id_2D_shift(idi1,1) = idi(idi1) + round(I/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_2D_shift(~idi1,1) = idi(~idi1) - floor(I/2);

idj1 = (idj <= floor(J/2));
id_2D_shift(idj1,2) = idj(idj1) + round(J/2); % Nfft/2 if even, (Nfft+1)/2 otherwise 
id_2D_shift(~idj1,2) = idj(~idj1) - floor(J/2);

% conversion to linear indices
id_shift = reshape((id_2D_shift(:,2) - 1)*I + id_2D_shift(:,1), size(id));

end

