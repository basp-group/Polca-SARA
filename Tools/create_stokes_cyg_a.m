 function [im] = create_stokes_cyg_a(im, stokes_P)
% Code to generate polarization images 
% x: Total intensity image
% p: Maximum fractional polarization
% stokes_P: Number of polarization images to be created (3: Stokes Q and U, 4: Stokes Q,U,V)

N = size(im{1});
per = 0.01;

im_max = max(im{1}(:));
im_th = zeros(N);

im_th(im{1}>0.1*im_max) = im{1}(im{1}>0.1*im_max);

im_true = im;

p_abs = sqrt(im{2}.^2+im{3}.^2);
p_err = p_abs - (im{1});
mask = zeros(N);
mask(p_err>0) = p_err(p_err>0);
p_count =  size(find(mask),1) ;

while p_count
    for i =2:stokes_P
%         im{i} = sign(im{i}).*(abs(im{i}) - p/2);
%           xx_p{i} = im{i} - (mask./2);
%           xx_n{i} = im{i} + (mask./2);
%           xx_p{i} = im{i}.*(1-per);
%           xx_n{i} = im{i}.*(1+per);
%           im{i}(im{i}>0) = xx_p{i}(im{i}>0);
%           im{i}(im{i}<0) = xx_n{i}(im{i}<0);
          im{i}(mask>0) = im{i}(mask>0).*(1-per);
        im{i}(abs(im{i})<10^(-3.4)) = 0; %0; %sign(im{i}(abs(im{i})<1e-5)).*1e-5; %1e-6
if i == 3
     im{i}(abs(im{i})<10^(-3.05)) = 0;
end
    end
    per = per+0.05;
    p_abs = sqrt(im{2}.^2+im{3}.^2);
    p_err = p_abs - (im{1});
    mask = zeros(N);
    mask(p_err>0) = p_err(p_err>0);
    p_count =  size(find(mask),1) ;
end

% param_im.im = im;

end