% Create 3C286 images from fits file

n1 = 100;
n2 = n1;

x1 = fitsread('3C286-I-16564.fits');
im{1} = imresize(x1,[n1,n2]);
% im{1}(im{1}<1e-3) = 1e-3;

x2 = fitsread('3C286-Q-16564.fits');
im{2} = imresize(x2,[n1,n2]);
% im{2}(im{2}<1e-4) = 1e-4;

x3 = fitsread('3C286-U-16564.fits');
im{3} = imresize(x3,[n1,n2]);
% im{3}(im{3}<1e-4) = 1e-4;