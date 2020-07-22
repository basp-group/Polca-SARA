% Create Stokes images for Cygnus A

n1 = 256; %256;
n2 = n1;
if n1 == 256
size_ker_h = [9,9]; %[9,9] ;
var_ker_h = 1; %1
else
    size_ker_h = [5,5];
    var_ker_h = 0.5;
end

h = fspecial('gaussian',size_ker_h,var_ker_h);

if hydra_ch == 1

x11 = fitsread('HYDRA-C-I.fits');
x1 = x11(258+50:769-50, 257+50:768-50);
im{1} = imresize(x1,[n1,n2]);

im_ = conv2(im{1},h,'same');
im{1} = (im_./max(im_(:))).*max(x1(:));
im{1}((im{1})<10^(-3.7)) = 0; 



 x21 = fitsread('HYDRA-C-Q.fits');
 x2 = x21(258+50:769-50, 257+50:768-50);
im{2} = imresize(x2,[n1,n2]);

im_ = conv2(im{2},h,'same');
im{2} = (im_./max(im_(:))).*max(x2(:));
im{2}(abs(im{2})<10^(-4.35)) = 0; 
im{2}(im{1} <= 0) = 0;


x31 = fitsread('HYDRA-C-U.fits');
x3 = x31(258+50:769-50, 257+50:768-50);
im{3} = imresize(x3,[n1,n2]);

 im_ = conv2(im{3},h,'same');
im{3} = (im_./max(im_(:))).*max(x3(:));
im{3}(abs(im{3})<10^(-4.5)) = 0; % 4.15
im{3}(im{1} <= 0) = 0;


im = create_stokes_hydra(im, stokes_P);

%%
% figure, imagesc(log10(abs(im_true{1}))), colorbar, axis image, axis off
% figure, imagesc(log10(abs(im_true{2}))), colorbar, axis image, axis off
% figure, imagesc(log10(abs(im_true{3}))), colorbar, axis image, axis off



%%
% Hydra Images for C band for A+B+C+D configuration for VLA

elseif hydra_ch == 2
x11 = fitsread('C-HYDRA-I.fits');
x1 = x11(544:1505, 544:1505);
im{1} = imresize(x1,[n1,n2]);

im_ = conv2(im{1},h,'same');
im{1} = (im_./max(im_(:))).*max(x1(:));
im{1}((im{1})<10^(-3.7)) = 0; 



 x21 = fitsread('C-HYDRA-Q.fits');
 x2 = x21(544:1505, 544:1505);
im{2} = imresize(x2,[n1,n2]);

im_ = conv2(im{2},h,'same');
im{2} = (im_./max(im_(:))).*max(x2(:));
im{2}(abs(im{2})<10^(-4.35)) = 0; 
im{2}(im{1} <= 0) = 0;


x31 = fitsread('C-HYDRA-U.fits');
x3 = x31(544:1505, 544:1505);
im{3} = imresize(x3,[n1,n2]);

 im_ = conv2(im{3},h,'same');
im{3} = (im_./max(im_(:))).*max(x3(:));
im{3}(abs(im{3})<10^(-4.15)) = 0; % 4.15
im{3}(im{1} <= 0) = 0;

im = create_stokes_hydra(im, stokes_P);




%%
% Hydra images for X band
elseif hydra_ch == 3
x11 = fitsread('HYDRA-X-I.fits');
x1 = x11(404:1645, 404:1645);
im{1} = imresize(x1,[n1,n2]);

im_ = conv2(im{1},h,'same');
im{1} = (im_./max(im_(:))).*max(x1(:));
im{1}((im{1})<10^(-3.7)) = 0; %3.1



 x21 = fitsread('HYDRA-X-Q.fits');
 x2 = x21(404:1645, 404:1645);
im{2} = imresize(x2,[n1,n2]);

im_ = conv2(im{2},h,'same');
im{2} = (im_./max(im_(:))).*max(x2(:));
im{2}(abs(im{2})<10^(-4.5)) = 0; %4.35
im{2}(im{1} <= 0) = 0;


x31 = fitsread('HYDRA-X-U.fits');
x3 = x31(404:1645, 404:1645);
im{3} = imresize(x3,[n1,n2]);

 im_ = conv2(im{3},h,'same');
im{3} = (im_./max(im_(:))).*max(x3(:));
im{3}(abs(im{3})<10^(-4.7)) = 0; % 4.46
im{3}(im{1} <= 0) = 0;

% 
% im{1} = 100*im{1};
im{2} = 2*im{2};
im{3} = 2*im{3};

im = create_stokes_hydra(im, stokes_P);

elseif hydra_ch == 5
x11 = fitsread('HYDRA-X-I.fits');
x1 = x11(404:1645, 404:1645);
im{1} = imresize(x1,[n1,n2]);

im_ = conv2(im{1},h,'same');
im{1} = (im_./max(im_(:))).*max(x1(:));
im{1}((im{1})<10^(-3.5)) = 0; %3.1



 x21 = fitsread('HYDRA-X-Q.fits');
 x2 = x21(404:1645, 404:1645);
im{2} = imresize(x2,[n1,n2]);

im_ = conv2(im{2},h,'same');
im{2} = (im_./max(im_(:))).*max(x2(:));
im{2}(abs(im{2})<10^(-4.5)) = 0; %4.35
im{2}(im{1} <= 0) = 0;


x31 = fitsread('HYDRA-X-U.fits');
x3 = x31(404:1645, 404:1645);
im{3} = imresize(x3,[n1,n2]);

 im_ = conv2(im{3},h,'same');
im{3} = (im_./max(im_(:))).*max(x3(:));
im{3}(abs(im{3})<10^(-4.6)) = 0; % 4.46
im{3}(im{1} <= 0) = 0;

% 
% im{1} = 100*im{1};
im{2} = 5*im{2};
im{3} = 5*im{3};

im = create_stokes_hydra(im, stokes_P);
end
