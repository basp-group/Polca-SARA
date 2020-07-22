function [im] = create_stokes_im(x,p, stokes_P)
% Code to generate polarization images 
% x: Total intensity image
% p: Maximum fractional polarization
%stokes_P: Number of polarization images to be created (3: Stokes Q and U, 4: Stokes Q,U,V)


per = 0.1;

if length(x) > 3
    N = size(x);
    im{1} = x;
    im_th = zeros(N);
    im_max = max(x(:));
    im_th(im{1}>0.1*im_max) = im{1}(im{1}>0.1*im_max);
for i =2:stokes_P
%     im{i}= zeros(N);
    im{i} = (p-0.1).*randn(N).*im_th;
end
else
    im = x;
    N = size(im{1});
end
im_true = im;

p_abs = sqrt(im{2}.^2+im{3}.^2);
p_err = p_abs - (p.*im{1});
mask = zeros(N);
mask(p_err>0) = p_err(p_err>0);
p_count =  size(find(mask),1) ;
t = 1;

while p_count
    for i =2:stokes_P
%           xx_p{i} = im{i}.*(1-per);
%           xx_n{i} = im{i}.*(1+per);
%           im{i}(im{i}>0) = xx_p{i}(im{i}>0);
%           im{i}(im{i}<0) = xx_n{i}(im{i}<0);
im{i}(p_err>0) = im{i}(p_err>0).*(1-per);
im{i}
    end
    per = per+0.05;
    p_abs = sqrt(im{2}.^2+im{3}.^2);
    p_err = p_abs - (p.*im{1});
    mask = zeros(N);
    mask(p_err>0) = p_err(p_err>0);
    p_count =  size(find(mask),1) ;
    s(t) = p_count;
    t = t+1;
    figure(100), semilogy(s)
    pause(0.1)
end

% param_im.im = im;

end