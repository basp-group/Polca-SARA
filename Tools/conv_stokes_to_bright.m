function [B] = conv_stokes_to_bright(St,L, Ni)
% Code to convert Stokes matrix into Brightness matrix
% St = [param_im.im_true{1}(:), param_im.im_true{2}(:), param_im.im_true{3}(:)];  % Stokes matrix 
 
 B1 = St*L; % Brightness matrix

 for i =1:4
    B{i} = reshape(B1(:,i), Ni(2), Ni(1)); 
 end
end
 